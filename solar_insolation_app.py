"""Streamlit app to explore Earth's orbit and daily insolation by latitude."""
import math
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib import font_manager

# ============================================
# 0. 날짜 → N일차
# ============================================
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def day_of_year(month: int, day: int) -> int:
    """Return zero-based day of year for the provided month/day."""
    safe_day = min(day, DAYS_IN_MONTH[month - 1])
    return sum(DAYS_IN_MONTH[: month - 1]) + safe_day - 1  # 0-based


def month_day_from_day_of_year(N: int) -> Tuple[int, int]:
    """Convert zero-based day-of-year to (month, day)."""

    remaining = N + 1  # convert to 1-based count for division into months
    for month_idx, days_in_month in enumerate(DAYS_IN_MONTH, start=1):
        if remaining > days_in_month:
            remaining -= days_in_month
        else:
            return month_idx, remaining
    # fallback: last day of the year
    return 12, 31


# ============================================
# 1. Mean anomaly
# ============================================
def mean_anomaly(N: int, M0_deg: float = -3.0) -> float:
    deg_per_day = 360.0 / 365.25
    M_deg = M0_deg + deg_per_day * N
    return math.radians(M_deg)


# ============================================
# 2. Eccentric anomaly
# ============================================
def eccentric_anomaly(M: float, e: float, n_iter: int = 6) -> float:
    E = M
    for _ in range(n_iter):
        f = E - e * math.sin(E) - M
        fprime = 1 - e * math.cos(E)
        E = E - f / fprime
    return E


# ============================================
# 3. True anomaly
# ============================================
def true_anomaly(E: float, e: float) -> float:
    num = math.sqrt(1 + e) * math.sin(E / 2)
    den = math.sqrt(1 - e) * math.cos(E / 2)
    v = 2 * math.atan2(num, den)
    return v


# ============================================
# 4. Declination (태양 적위)
# ============================================
def solar_declination(v: float, omega_deg: float, epsilon_deg: float) -> Tuple[float, float]:
    # 지구의 진근점 기준 경도 v에 근일점 경도(omega)를 더한 뒤 태양-지구 시차 180°를
    # 반영해 실제 태양 황경(lam)을 얻는다. 이렇게 하면 북반구 하지가 원일점,
    # 동지가 근일점에 일치한다.
    lam = v + math.radians(omega_deg) + math.pi
    eps = math.radians(epsilon_deg)
    delta = math.asin(math.sin(eps) * math.sin(lam))
    return delta, lam


# ============================================
# 5. 위도별 일사량
# ============================================
S0 = 1361  # solar constant


def ensure_korean_font() -> None:
    """Ensure matplotlib uses a font that can render Korean labels."""

    preferred_fonts = [
        "NanumGothic",
        "Noto Sans CJK KR",
        "Malgun Gothic",
        "AppleGothic",
    ]

    for font_name in preferred_fonts:
        try:
            font_manager.findfont(font_name, fallback_to_default=False)
        except Exception:
            continue
        else:
            plt.rcParams["font.family"] = font_name
            break
    else:
        font_url = (
            "https://github.com/google/fonts/raw/main/ofl/"
            "nanumgothic/NanumGothic-Regular.ttf"
        )
        cache_path = Path(tempfile.gettempdir()) / "NanumGothic-Regular.ttf"
        try:
            if not cache_path.exists():
                urllib.request.urlretrieve(font_url, cache_path)
            font_manager.fontManager.addfont(cache_path)
            plt.rcParams["font.family"] = "NanumGothic"
        except Exception:
            plt.rcParams["font.family"] = "DejaVu Sans"

    plt.rcParams["axes.unicode_minus"] = False


def daily_insolation(phi_rad: np.ndarray, delta: float, e: float, lam: float, omega_deg: float) -> np.ndarray:
    rfac = ((1 + e * math.cos(lam - math.radians(omega_deg))) / (1 - e * e)) ** 2

    cosH0 = -np.tan(phi_rad) * np.tan(delta)
    cosH0 = np.clip(cosH0, -1, 1)
    H0 = np.arccos(cosH0)

    Q = (S0 / math.pi) * rfac * (
        H0 * np.sin(phi_rad) * np.sin(delta)
        + np.cos(phi_rad) * np.cos(delta) * np.sin(H0)
    )
    return Q


# ============================================
# 6. 남중고도
# ============================================
def solar_noon_altitude(phi_rad: float, delta: float) -> float:
    alpha = math.degrees(
        math.asin(
            math.sin(phi_rad) * math.sin(delta)
            + math.cos(phi_rad) * math.cos(delta)
        )
    )
    return alpha


# ============================================
# 7. 궤도 그림
# ============================================
def draw_orbit(e: float, omega_deg: float, E_now: float, epsilon_deg: float):
    # 궤도 시각화에서만 이심률을 과장해 학생들이 타원 형태를 더 쉽게 구분하도록 한다.
    e_vis = min(e * 10, 0.9)  # 시각화용 이심률 (계산은 실제 e 사용)

    a = 1.0
    b = a * math.sqrt(1 - e_vis * e_vis)
    omega = math.radians(omega_deg)
    eps = math.radians(epsilon_deg)

    E_all = np.linspace(0, 2 * np.pi, 500)
    x = a * (np.cos(E_all) - e_vis)
    y = b * np.sin(E_all)

    # 회전
    xR = x * np.cos(omega) - y * np.sin(omega)
    yR = x * np.sin(omega) + y * np.cos(omega)

    # 현재 지구 위치
    xE = a * (math.cos(E_now) - e_vis)
    yE = b * math.sin(E_now)
    xE_R = xE * np.cos(omega) - yE * np.sin(omega)
    yE_R = xE * np.sin(omega) + yE * np.cos(omega)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(xR, yR, label="궤도")
    ax.scatter(0, 0, s=200, color="yellow", label="태양")

    # Earth
    ax.scatter(xE_R, yE_R, color="blue", s=100, label="지구")

    # 자전축 (2D)
    L = 0.3
    dx = L * math.sin(eps)
    dy = L * math.cos(eps)
    ax.plot([xE_R - dx / 2, xE_R + dx / 2], [yE_R - dy / 2, yE_R + dy / 2], color="black", linewidth=2)

    ax.set_aspect("equal")
    R = 1 + e + 0.5
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.grid()
    ax.legend()
    ax.set_title(
        f"지구 공전 궤도 (시각화용 e={e_vis:.3f}, 실제 e={e:.4f})"
    )
    return fig


# ============================================
# STREAMLIT APP (학생용 UI)
# ============================================
st.set_page_config(layout="wide")
st.title("🌍 고등학생용 지구 공전·일사량 시뮬레이터")
ensure_korean_font()
st.caption(
    "달력 날짜와 공전 매개변수를 조정하며 일사량 변화를 한눈에 확인하세요. "
    "왼쪽 사이드바에서 값을 바꾼 뒤, 아래 애니메이션으로 흐름을 살펴볼 수 있습니다."
)

# --- session_state ---
if "animate" not in st.session_state:
    st.session_state.animate = False
if "N" not in st.session_state:
    st.session_state.N = 80

# --------------------------------------------
# 입력 UI
# --------------------------------------------
with st.sidebar:
    st.header("입력값")

    st.subheader("날짜 선택")
    date_mode = st.radio("날짜 입력 방식", ("월·일로 입력", "N일차 슬라이더"), index=0)

    if date_mode == "월·일로 입력":
        month = int(st.number_input("월", 1, 12, 3))
        max_day_for_month = DAYS_IN_MONTH[month - 1]
        day = int(
            st.number_input(
                "일", 1, max_day_for_month, min(21, max_day_for_month)
            )
        )

        st.caption(f"현재 달에서 계산되는 최대 날짜는 {max_day_for_month}일입니다.")
        N_from_date = min(day_of_year(month, day), 364)
        N_slider = None
    else:
        month, day = None, None
        N_slider = st.slider("날짜(N일차)", 0, 364, st.session_state.N)
        N_from_date = None

    st.subheader("공전 매개변수")
    e = st.slider("이심률 e", 0.0, 0.1, 0.0167, 0.0001)
    omega_deg = st.slider("세차(ω)", 0.0, 360.0, 102.0)
    epsilon_deg = st.slider("축 경사(ε)", 0.0, 40.0, 23.4)

    st.subheader("위치·시간")
    phi_deg = st.slider("위도", -90.0, 90.0, 37.0)

    st.divider()
    anim_speed = st.slider("애니메이션 속도 (ms)", 1, 250, 30, 1)

if date_mode == "월·일로 입력" and st.session_state.animate:
    st.session_state.animate = False

# --------------------------------------------
# 레이아웃
# --------------------------------------------
colL, colR = st.columns([1.15, 1])

with colL:
    if not st.session_state.animate and date_mode == "N일차 슬라이더" and N_slider is not None:
        st.session_state.N = N_slider

    if date_mode == "월·일로 입력":
        active_N = N_from_date
    else:
        active_N = st.session_state.N

    M = mean_anomaly(active_N)
    E_val = eccentric_anomaly(M, e)
    v = true_anomaly(E_val, e)
    delta, lam = solar_declination(v, omega_deg, epsilon_deg)

    fig_orbit = draw_orbit(e, omega_deg, E_val, epsilon_deg)
    st.pyplot(fig_orbit, width="stretch")

with colR:
    st.subheader("🌞 선택 날짜와 태양 위치")
    if date_mode == "월·일로 입력":
        st.markdown(f"**입력한 날짜:** {month}월 {day}일")
    else:
        derived_month, derived_day = month_day_from_day_of_year(active_N)
        st.markdown(
            f"**슬라이더 N일차:** {active_N}일차 · **달력 환산:** {derived_month}월 {derived_day}일"
        )

    st.markdown(
        f"**위도:** {phi_deg:.1f}° · **태양 적위:** {math.degrees(delta):.2f}°"
    )
    phi_rad = math.radians(phi_deg)
    alpha_noon = solar_noon_altitude(phi_rad, delta)

    st.subheader("📈 위도별 하루 태양 에너지량")
    phi_list = np.linspace(-90, 90, 181)
    phi_rad_all = np.radians(phi_list)
    Q = daily_insolation(phi_rad_all, delta, e, lam, omega_deg)

    figQ, axQ = plt.subplots(figsize=(6, 4))
    axQ.plot(phi_list, Q, color="#1f77b4", linewidth=2.2)
    axQ.fill_between(phi_list, Q, color="#1f77b4", alpha=0.08)
    axQ.axvline(phi_deg, color="crimson", linestyle="--", linewidth=1.5, label="선택 위도")
    axQ.legend()
    axQ.set_xlabel("위도 (deg)")
    axQ.set_ylabel("태양 에너지량 (W/m²)")
    axQ.spines["top"].set_visible(False)
    axQ.spines["right"].set_visible(False)
    axQ.grid(alpha=0.3)
    st.pyplot(figQ, width="stretch")

    st.subheader("🌅 남중고도")
    st.metric("정오 고도", f"{alpha_noon:.2f}°")

    st.divider()
    st.markdown(
        """
        - 그래프에 면적 음영과 범례를 추가해 위도 변화에 따른 일사량의 흐름을 쉽게 읽을 수 있습니다.
        """
    )

# --------------------------------------------
# 애니메이션 컨트롤
# --------------------------------------------
st.subheader("⏯ 날짜 자동 변화 애니메이션")

if date_mode == "월·일로 입력":
    st.info("월·일 입력 모드에서는 N일차가 자동 계산됩니다. 슬라이더 모드에서 애니메이션을 사용할 수 있습니다.")
    c1, c2 = st.columns(2)
    c1.button("▶ Start", disabled=True)
    c2.button("⏸ Pause", disabled=True)
else:
    c1, c2, c3 = st.columns(3)
    if c1.button("▶ Start"):
        st.session_state.animate = True
    if c2.button("⏸ Pause"):
        st.session_state.animate = False
    if c3.button("↺ 0일차로 리셋"):
        st.session_state.N = 0

if st.session_state.animate:
    st.session_state.N = (st.session_state.N + 1) % 365
    time.sleep(anim_speed / 1000.0)
    st.experimental_rerun()
