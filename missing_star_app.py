"""Streamlit app to generate zenith star charts with randomly removed stars."""
from __future__ import annotations

import math
import os
import random
from datetime import datetime
from typing import Iterable, Sequence
from zoneinfo import ZoneInfo

import streamlit as st
from starplot import Observer, Star, ZenithPlot, _

# 최대 표시 등급
MAX_PLOT_MAG = 4.0
# 이미지 출력 경로
CHART_DIR = "charts"
PROBLEM_PATH = os.path.join(CHART_DIR, "problem.png")
ANSWER_PATH = os.path.join(CHART_DIR, "answer.png")
DEFAULT_TZ = ZoneInfo("Asia/Seoul")


# =====================================================
# 유틸 함수
# =====================================================

def calc_alt_deg(star: Star, obs: Observer) -> float:
    """Calculate altitude (deg) for a star given an observer."""
    lat_rad = math.radians(obs.lat)
    dec_rad = math.radians(star.dec)

    lst_deg = obs.lst
    ha_deg = (lst_deg - star.ra) % 360
    ha_rad = math.radians(ha_deg)

    sin_alt = (
        math.sin(lat_rad) * math.sin(dec_rad)
        + math.cos(lat_rad) * math.cos(dec_rad) * math.cos(ha_rad)
    )
    sin_alt = max(-1.0, min(1.0, sin_alt))
    return math.degrees(math.asin(sin_alt))


def stars_above_horizon(stars: Sequence[Star], obs: Observer) -> list[Star]:
    return [s for s in stars if calc_alt_deg(s, obs) > 0]


def pick_missing_stars(candidates: Sequence[Star], k: int) -> list[Star]:
    return random.sample(list(candidates), k)


def ensure_chart_dir() -> None:
    os.makedirs(CHART_DIR, exist_ok=True)


# =====================================================
# Streamlit UI
# =====================================================


def render_form() -> dict:
    col1, col2 = st.columns(2)
    with col1:
        date_input = st.date_input("날짜 선택", datetime.now(tz=DEFAULT_TZ).date())
        time_input = st.time_input(
            "시간 선택", value=datetime.now(tz=DEFAULT_TZ).time()
        )
    with col2:
        lat = st.number_input("위도 입력", value=37.5665, format="%.6f")
        lon = st.number_input("경도 입력", value=126.9780, format="%.6f")

    n = st.number_input(
        "삭제 후보 최대 등급 n",
        value=3.0,
        min_value=0.0,
        max_value=MAX_PLOT_MAG,
        step=0.1,
    )
    k = st.number_input("삭제할 별 수 k", value=3, min_value=1, step=1)
    return {
        "date": date_input,
        "time": time_input,
        "lat": lat,
        "lon": lon,
        "n": n,
        "k": int(k),
    }


def build_observer(date_val: datetime.date, time_val, lat: float, lon: float) -> Observer:
    dt = datetime.combine(date_val, time_val).replace(tzinfo=DEFAULT_TZ)
    return Observer(dt=dt, lat=lat, lon=lon)


def select_candidates(n: float) -> Iterable[Star]:
    return Star.find(where=[_.magnitude <= MAX_PLOT_MAG, _.magnitude <= n, _.hip.notnull()])


def make_problem_plot(observer: Observer, missing_hip_ids: set[int]) -> ZenithPlot:
    plot = ZenithPlot(observer=observer, resolution=3000, scale=0.9)
    if missing_hip_ids:
        hip_list = ",".join(str(h) for h in missing_hip_ids)
        sql_filter = (
            "select * from _ "
            f"where magnitude <= {MAX_PLOT_MAG} "
            f"and (hip is null or hip not in ({hip_list}))"
        )
    else:
        sql_filter = f"select * from _ where magnitude <= {MAX_PLOT_MAG}"
    plot.stars(sql=sql_filter, where_labels=[False])
    plot.horizon()
    return plot


def make_answer_plot(observer: Observer, missing_hip_ids: set[int]) -> ZenithPlot:
    plot = ZenithPlot(observer=observer, resolution=3000, scale=0.9)
    plot.constellations()
    plot.stars(where=[_.magnitude <= MAX_PLOT_MAG], where_labels=[False])
    if missing_hip_ids:
        plot.stars(
            where=[_.hip.isin(list(missing_hip_ids))],
            where_labels=[False],
            style__marker__color="red",
            style__marker__size=18,
        )
    plot.horizon()
    return plot


def export_plot(plot: ZenithPlot, path: str) -> None:
    plot.export(path, transparent=True)


def render_results(problem_path: str, answer_path: str, missing_stars: Sequence[Star]) -> None:
    st.success("성도 생성 완료!")
    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("문제 성도")
        st.image(problem_path)
    with col_b:
        st.subheader("정답 성도")
        st.image(answer_path)
    st.subheader("삭제된 별 목록 (HIP / 등급)")
    st.write([f"HIP {s.hip} | mag={s.magnitude:.2f}" for s in missing_stars])


def main() -> None:
    st.set_page_config(page_title="Missing Star Generator", layout="wide")
    st.title("⭐ 미싱 스타 성도 생성기 (Streamlit)")
    st.write("날짜/시간, 위치, 밝기 등급을 선택하면 자동으로 문제/정답 성도를 만들어줍니다.")

    params = render_form()
    run_btn = st.button("👉 성도 생성하기")

    if not run_btn:
        return

    ensure_chart_dir()
    observer = build_observer(params["date"], params["time"], params["lat"], params["lon"])

    candidates = list(select_candidates(params["n"]))
    candidates = stars_above_horizon(candidates, observer)
    k = params["k"]
    if len(candidates) < k:
        st.error(f"지평선 위의 삭제 후보 별이 {len(candidates)}개인데 k={k}개를 요청했습니다.")
        return

    missing_stars = pick_missing_stars(candidates, k)
    missing_hip_ids = {s.hip for s in missing_stars}

    problem_plot = make_problem_plot(observer, missing_hip_ids)
    answer_plot = make_answer_plot(observer, missing_hip_ids)

    export_plot(problem_plot, PROBLEM_PATH)
    export_plot(answer_plot, ANSWER_PATH)

    render_results(PROBLEM_PATH, ANSWER_PATH, missing_stars)


if __name__ == "__main__":
    main()
