"""Streamlit app to explore Earth's orbit and daily insolation by latitude."""

import math
import time
import urllib.request
from pathlib import Path
from typing import Tuple
import tempfile
import os

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib import font_manager

# ============================================
# 0. 날짜 처리
# ============================================
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
YEAR_DAYS = 365.2422


def trigger_rerun():
    """Use Streamlit's stable rerun API with a safe fallback for older versions."""
    if hasattr(st, "rerun"):
        st.rerun()
    else:  # pragma: no cover - legacy Streamlit
        st.experimental_rerun()


def ensure_korean_font():
    """matplotlib에서 한글이 깨지지 않도록 강제로 설정."""
    font_path = "/tmp/NanumGothic.ttf"
    url = (
        "https://github.com/google/fonts/raw/main/ofl/"
        "nanumgothic/NanumGothic-Regular.ttf"
    )

    if not os.path.exists(font_path):
        urllib.request.urlretrieve(url, font_path)

    font_manager.fontManager.addfont(font_path)
    plt.rcParams["font.family"] = "NanumGothic"
    plt.rcParams["axes.unicode_minus"] = False


def day_of_year(month: int, day: int) -> int:
    safe_day = min(day, DAYS_IN_MONTH[month - 1])
    return sum(DAYS_IN_MONTH[: month - 1]) + safe_day - 1  # 0-based


def month_day_from_day_of_year(N: int) -> Tuple[int, int]:
    """Zero-based N → (month, day)"""
    remaining = N + 1  # convert to 1-based day
    for month_idx, days_in_month in enumerate(DAYS_IN_MONTH, start=1):
        if remaining > days_in_month:
            remaining -= days_in_month
        else:
            return month_idx, remaining
    return 12, 31


# ============================================
# 1. Orbital mechanics
# ============================================
def eccentric_from_true(v: float, e: float) -> float:
    """Return eccentric anomaly for a given true anomaly."""

    factor = math.sqrt((1 - e) / (1 + e))
    return 2 * math.atan2(factor * math.tan(v / 2), 1)


def eccentric_anomaly(M: float, e: float, n_iter: int = 6) -> float:
    """Solve Kepler's equation M = E - e sin(E) for eccentric anomaly E."""

    E = M
    for _ in range(n_iter):
        f = E - e * math.sin(E) - M
        fprime = 1 - e * math.cos(E)
        E = E - f / fprime
    return E


def true_anomaly_from_eccentric(E: float, e: float) -> float:
    num = math.sqrt(1 + e) * math.sin(E / 2)
    den = math.sqrt(1 - e) * math.cos(E / 2)
    return 2 * math.atan2(num, den)


def solar_declination(lam: float, epsilon_deg: float) -> float:
    """Return solar declination from ecliptic longitude."""

    eps = math.radians(epsilon_deg)
    return math.asin(math.sin(eps) * math.sin(lam))


# ============================================
# 2. Insolation
# ============================================
S0 = 1361


def daily_insolation(phi_rad: np.ndarray, delta: float, e: float, lam: float, omega_deg: float) -> np.ndarray:
    omega_rad = math.radians(omega_deg)
    rfac = ((1 + e * math.cos(lam - omega_rad)) / (1 - e * e)) ** 2

    cosH0 = -np.tan(phi_rad) * np.tan(delta)
    cosH0 = np.clip(cosH0, -1, 1)
    H0 = np.arccos(cosH0)

    return (S0 / math.pi) * rfac * (
        H0 * np.sin(phi_rad) * np.sin(delta)
        + np.cos(phi_rad) * np.cos(delta) * np.sin(H0)
    )


def solar_noon_altitude(phi_rad: float, delta: float) -> float:
    return math.degrees(
        math.asin(
            math.sin(phi_rad) * math.sin(delta) +
            math.cos(phi_rad) * math.cos(delta)
        )
    )


def sun_path_curve(phi_rad: float, delta: float) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return hour offsets (radians), altitude (radians), and noon altitude."""
    cosH0 = -math.tan(phi_rad) * math.tan(delta)
    cosH0 = max(-1.0, min(1.0, cosH0))
    H0 = math.acos(cosH0)

    H_range = np.linspace(-H0, H0, 200)
    alt = np.arcsin(
        np.sin(phi_rad) * np.sin(delta) +
        np.cos(phi_rad) * np.cos(delta) * np.cos(H_range)
    )

    noon_alt = float(np.max(alt))
    return H_range, alt, noon_alt


# ============================================
# 3. Orbit visualization
# ============================================
def draw_orbit(e: float, omega_deg: float, E_now: float, epsilon_deg: float):
    """공전 위상·세차·자전축을 물리 순서대로 적용한 궤도 시각화."""

    # 시각적으로 보기 좋게 이심률 강조
    e_vis = min(e * 10, 0.9)

    # 기본 파라미터
    a = 1.0
    b = a * math.sqrt(1 - e_vis * e_vis)
    omega_rad = math.radians(omega_deg)
    eps_rad = math.radians(epsilon_deg)
    view_tilt = math.radians(15)

    # 회전 행렬
    def R_z(theta: float) -> np.ndarray:
        return np.array(
            [
                [math.cos(theta), -math.sin(theta), 0],
                [math.sin(theta), math.cos(theta), 0],
                [0, 0, 1],
            ]
        )

    def R_x(theta: float) -> np.ndarray:
        return np.array(
            [
                [1, 0, 0],
                [0, math.cos(theta), -math.sin(theta)],
                [0, math.sin(theta), math.cos(theta)],
            ]
        )

    view_rot = R_x(-view_tilt)

    # 공전 궤도 좌표 (perihelion이 왼쪽)
    E_all = np.linspace(0, 2 * np.pi, 500)
    x_base = -a * (np.cos(E_all) - e_vis)
    y_base = b * np.sin(E_all)
    orbit_base = np.vstack([x_base, y_base, np.zeros_like(x_base)])
    orbit_rot = R_z(omega_rad) @ orbit_base
    orbit_view = view_rot @ orbit_rot
    xR, yR = orbit_view[0], orbit_view[1]

    # 현재 지구 위치
    xE_base = -a * (math.cos(E_now) - e_vis)
    yE_base = b * math.sin(E_now)
    pos_base = np.array([xE_base, yE_base, 0.0])
    pos_rot = R_z(omega_rad) @ pos_base
    pos_view = view_rot @ pos_rot
    xE_R, yE_R = pos_view[0], pos_view[1]

    # 근일점 / 원일점 (세차에 따라 회전)
    peri_base = np.array([-a * (1 - e_vis), 0.0, 0.0])
    ap_base = np.array([a * (1 + e_vis), 0.0, 0.0])
    peri_rot = R_z(omega_rad) @ peri_base
    ap_rot = R_z(omega_rad) @ ap_base
    peri_view = view_rot @ peri_rot
    ap_view = view_rot @ ap_rot

    # 자전축: 우주 공간에서 방향 고정 (세차·경사만 적용)
    base_axis = np.array([0.0, 0.0, 1.0])
    axis_3d = R_z(omega_rad) @ R_x(-eps_rad) @ base_axis
    axis_view = view_rot @ axis_3d
    axis_proj = axis_view[:2] / np.linalg.norm(axis_view[:2])

    L = 0.35
    xA1, yA1 = xE_R - axis_proj[0] * L, yE_R - axis_proj[1] * L
    xA2, yA2 = xE_R + axis_proj[0] * L, yE_R + axis_proj[1] * L

    # 그림
    fig, ax = plt.subplots(figsize=(2.4, 2.4), facecolor="#0a0f1c")
    ax.set_facecolor("#0a0f1c")

    ax.plot(xR, yR, linestyle="--", color="#4db7ff", linewidth=1.6)
    ax.scatter(0, 0, s=140, color="#0d5c84", edgecolors="white", linewidths=1.5)
    ax.text(0, 0, "☀", color="#ffef9f", fontsize=18, ha="center", va="center", weight="bold")
    ax.scatter(xE_R, yE_R, s=70, color="#78ffba", edgecolors="#0a0f1c", linewidths=1.2)
    ax.plot([xA1, xA2], [yA1, yA2], color="white", linewidth=2)
    ax.scatter(peri_view[0], peri_view[1], s=70, color="#78ffba", alpha=0.9)
    ax.scatter(ap_view[0], ap_view[1], s=70, color="#78ffba", alpha=0.9)

    ax.set_aspect("equal")
    R = 1 + e_vis + 0.55
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.tick_params(colors="#0a0f1c", labelsize=6)  # hide coordinates
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#1c2f46")
    ax.grid(color="#1c2f46", linestyle="--", linewidth=0.6)
    fig.tight_layout(pad=0.6)

    return fig


def draw_sun_path(phi_deg: float, delta: float, epsilon_deg: float):
    """2D 천구도 스타일로 태양 고도 변화를 표시."""

    phi_rad = math.radians(phi_deg)
    eps_rad = math.radians(epsilon_deg)
    declinations = [eps_rad, 0.0, -eps_rad]
    labels = ["하지", "춘/추분", "동지"]
    colors = ["#f5c542", "#6fb3ff", "#9ad7a8"]

    fig, ax = plt.subplots(figsize=(2.8, 2.6))

    ax.axhspan(0, 90, color="#eef2ff", alpha=0.8)
    ax.axhline(0, color="#9ca3af", linewidth=1.0)
    ax.text(12, -7, "지평선", ha="center", va="top", fontsize=8, color="#4b5563")

    for dec, label, color in zip(declinations, labels, colors):
        H_range, alt, noon_alt = sun_path_curve(phi_rad, dec)
        hours = 12 + (H_range / math.pi) * 12
        alt_deg = np.degrees(alt)
        ax.plot(hours, alt_deg, color=color, linewidth=1.4, label=f"{label} (δ={math.degrees(dec):.1f}°)")
        ax.scatter(12, math.degrees(noon_alt), color=color, s=20, zorder=3)

    H_sel, alt_sel, noon_alt_sel = sun_path_curve(phi_rad, delta)
    hours_sel = 12 + (H_sel / math.pi) * 12
    alt_sel_deg = np.degrees(alt_sel)
    ax.plot(hours_sel, alt_sel_deg, color="#ff6b6b", linewidth=2.0, label="선택 날짜")
    ax.scatter(12, math.degrees(noon_alt_sel), color="#ff6b6b", s=26, zorder=4)

    ax.set_xlim(0, 24)
    ax.set_ylim(-10, 95)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_yticks(range(0, 91, 30))
    ax.set_xlabel("태양시 (h)", fontsize=9)
    ax.set_ylabel("고도 (°)", fontsize=9)
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    ax.set_title("하늘에서 본 태양 경로", fontsize=11)
    ax.grid(alpha=0.3)
    fig.tight_layout(pad=0.4)
    return fig


# ============================================
# STREAMLIT APP (학생용)
# ============================================
st.set_page_config(layout="wide")
ensure_korean_font()

st.title("밀란코비치 주기에 따른 기후 변화")

# 상단 여백 축소 및 헤더 간격 조정
st.markdown(
    """
    <style>
    .block-container { padding-top: 1.3rem; padding-bottom: 1.3rem; }
    h1 { margin-bottom: 0.4rem; }
    body { background-color: white; }
    </style>
    """,
    unsafe_allow_html=True,
)

INIT_MONTH = 3
INIT_DAY = 21
INIT_E = 0.0167
INIT_PRECESSION_YEAR = 0
INIT_OMEGA = (INIT_PRECESSION_YEAR / 26000) * 360
INIT_EPS = 23.44
INIT_PHI = 37.0
INIT_N = day_of_year(INIT_MONTH, INIT_DAY)
INIT_SPEED = 30
EQUINOX_N = day_of_year(3, 20)

if "animate" not in st.session_state:
    st.session_state.animate = False
if "N" not in st.session_state:
    st.session_state.N = INIT_N
if "month" not in st.session_state:
    st.session_state.month = INIT_MONTH
if "day" not in st.session_state:
    st.session_state.day = INIT_DAY
if "anim_speed" not in st.session_state:
    st.session_state.anim_speed = INIT_SPEED
if "phi_deg" not in st.session_state:
    st.session_state.phi_deg = INIT_PHI
if "e" not in st.session_state:
    st.session_state.e = INIT_E
if "precession_year" not in st.session_state:
    st.session_state.precession_year = INIT_PRECESSION_YEAR
if "omega_deg" not in st.session_state:
    st.session_state.omega_deg = INIT_OMEGA
if "epsilon_deg" not in st.session_state:
    st.session_state.epsilon_deg = INIT_EPS


def reset_state():
    """Restore all interactive values to the app's initial defaults."""
    st.session_state.month = INIT_MONTH
    st.session_state.day = INIT_DAY
    st.session_state.N = INIT_N
    st.session_state.animate = False
    st.session_state.e = INIT_E
    st.session_state.precession_year = INIT_PRECESSION_YEAR
    st.session_state.omega_deg = INIT_OMEGA
    st.session_state.epsilon_deg = INIT_EPS
    st.session_state.phi_deg = INIT_PHI
    st.session_state.anim_speed = INIT_SPEED

# --------------------------------------------
# 입력 UI (사이드바)
# --------------------------------------------
with st.sidebar:
    st.subheader("날짜 선택")

    # 애니메이션 중에는 month/day 자동 갱신
    if st.session_state.animate:
        m, d = month_day_from_day_of_year(st.session_state.N)
        st.session_state.month = m
        st.session_state.day = d

    date_cols = st.columns(2)
    month = int(
        date_cols[0].selectbox("월", list(range(1, 13)), index=st.session_state.month - 1)
    )
    max_day = DAYS_IN_MONTH[month - 1]
    day = int(
        date_cols[1].selectbox(
            "일",
            list(range(1, max_day + 1)),
            index=min(st.session_state.day, max_day) - 1,
        )
    )
    st.session_state.month = month
    st.session_state.day = day

    shortcuts = st.columns(4)
    if shortcuts[0].button("춘분"):
        st.session_state.month = 3
        st.session_state.day = 20
        st.session_state.animate = False
        trigger_rerun()

    if shortcuts[1].button("하지"):
        st.session_state.month = 6
        st.session_state.day = 21
        st.session_state.animate = False
        trigger_rerun()

    if shortcuts[2].button("추분"):
        st.session_state.month = 9
        st.session_state.day = 22
        st.session_state.animate = False
        trigger_rerun()

    if shortcuts[3].button("동지"):
        st.session_state.month = 12
        st.session_state.day = 21
        st.session_state.animate = False
        trigger_rerun()

    st.subheader("관측자 위도")
    lat_cols = st.columns([1, 2.4, 1])
    with lat_cols[0]:
        if st.button("−", help="위도 1° 감소"):
            st.session_state.phi_deg = max(-90.0, st.session_state.phi_deg - 1)
    with lat_cols[1]:
        phi_deg = st.number_input(
            "위도 (°)",
            -90.0,
            90.0,
            float(st.session_state.phi_deg),
            step=0.5,
            format="%.1f",
        )
        st.session_state.phi_deg = phi_deg
    with lat_cols[2]:
        if st.button("+", help="위도 1° 증가"):
            st.session_state.phi_deg = min(90.0, st.session_state.phi_deg + 1)

    st.subheader("밀란코비치 변수")
    e = st.slider("이심률 e", 0.0, 0.1, 0.0167, 0.0001, key="e")
    precession_year = st.slider(
        "세차 단계 (년)",
        0,
        26000,
        st.session_state.precession_year,
        key="precession_year",
        help=(
            "지구 세차운동은 약 26000년 주기로 반시계 방향으로 진행하며,\n"
            "계절의 위치(근일점/원일점과의 상대 위치)가 서서히 변합니다."
        ),
    )
    omega_deg = (precession_year / 26000) * 360
    st.session_state.omega_deg = omega_deg
    epsilon_deg = st.slider("축 경사(ε)", 0.0, 40.0, 23.44, key="epsilon_deg")

    st.subheader("애니메이션 속도")
    st.session_state.anim_speed = st.slider("속도(ms)", 1, 200, st.session_state.anim_speed)

# --------------------------------------------
# 메인 패널
# --------------------------------------------
# 날짜 → N 변환
active_N = (
    day_of_year(st.session_state.month, st.session_state.day)
    if not st.session_state.animate
    else st.session_state.N
)

omega_rad = math.radians(omega_deg)

# 날짜 → 평균황경(L) → 평균근점이각(M=L-ω) → 편심근점이각 → 진근점이각 → 황경
mean_longitude = 2 * math.pi * (active_N / YEAR_DAYS)
M = (mean_longitude - omega_rad) % (2 * math.pi)
E_val = eccentric_anomaly(M, e)
v = true_anomaly_from_eccentric(E_val, e)
lam = (v + omega_rad) % (2 * math.pi)
delta = solar_declination(lam, epsilon_deg)

phi_list = np.linspace(-90, 90, 181)
phi_rad_all = np.radians(phi_list)
Q = daily_insolation(phi_rad_all, delta, e, lam, omega_deg)

phi_rad = math.radians(phi_deg)
alpha_noon = solar_noon_altitude(phi_rad, delta)
Q_at_lat = float(daily_insolation(np.array([phi_rad]), delta, e, lam, omega_deg)[0])

cosH0 = -math.tan(phi_rad) * math.tan(delta)
cosH0 = max(-1.0, min(1.0, cosH0))
H0 = math.acos(cosH0)
daylight_hours = 24 * H0 / math.pi
hours = int(daylight_hours)
minutes = round((daylight_hours - hours) * 60)
if minutes == 60:
    hours += 1
    minutes = 0

# 단순한 평균 기온 추정 (위도와 계절 위상 기반의 학습용 모델)
season_phase = 2 * math.pi * (active_N - 80) / 365.0
base_temp = 15 - (abs(phi_deg) / 90.0) * 30
seasonal_amp = 10 * math.sqrt(max(math.cos(math.radians(phi_deg)), 0))
avg_temp = base_temp + seasonal_amp * math.sin(season_phase)

info_table_html = f"""
<style>
.info-table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 14px;
  table-layout: fixed;
}}
.info-table th {{
  background: #f3f4f6;
  color: #111827;
  padding: 8px 10px;
  font-weight: 800;
  text-align: left;
  border-bottom: 1px solid #e5e7eb;
}}
.info-table td {{
  padding: 8px 10px;
  border-bottom: 1px solid #e5e7eb;
  color: #111827;
  font-weight: 700;
  text-align: center;
}}
.info-table td.value {{ text-align: center; color: #0f172a; }}
</style>
<table class="info-table">
  <tr>
    <th>입력 날짜</th>
    <th>위도</th>
    <th>태양 적위</th>
    <th>정오 고도</th>
    <th>일사량</th>
    <th>낮 길이</th>
    <th>평균 기온</th>
  </tr>
  <tr>
    <td class="value">{st.session_state.month}월 {st.session_state.day}일</td>
    <td class="value">{phi_deg:.1f}°</td>
    <td class="value">{math.degrees(delta):.2f}°</td>
    <td class="value">{alpha_noon:.2f}°</td>
    <td class="value">{Q_at_lat:.0f} W/m²</td>
    <td class="value">{hours}시간 {minutes:02d}분</td>
    <td class="value">{avg_temp:.1f}°C</td>
  </tr>
</table>
"""

top_col_orbit, top_col_chart, top_col_sky = st.columns([1, 1, 1])

with top_col_orbit:
    st.subheader("🛰️ 지구 공전 궤도")
    fig_orbit = draw_orbit(e, omega_deg, E_val, epsilon_deg)
    st.pyplot(fig_orbit)

    ctrl_cols = st.columns([1, 1, 1])

    with ctrl_cols[0]:
        if st.button("▶", help="재생"):
            st.session_state.N = day_of_year(st.session_state.month, st.session_state.day)
            st.session_state.animate = True

    with ctrl_cols[1]:
        if st.button("⏸", help="일시정지"):
            st.session_state.animate = False

    with ctrl_cols[2]:
        if st.button("↺", help="처음 상태로 초기화"):
            reset_state()
            trigger_rerun()

with top_col_chart:
    st.subheader("📈 위도별 하루 태양 에너지량 (W/m²)")

    figQ, axQ = plt.subplots(figsize=(2.8, 2.8))
    axQ.plot(phi_list, Q, linewidth=1.6)
    axQ.axvline(phi_deg, color="red", linestyle="--", linewidth=1.2)
    axQ.grid(alpha=0.3)
    axQ.set_xlabel("위도", fontsize=9)
    axQ.set_ylabel("일사량(W/m²)", fontsize=9)

    xticks = np.arange(-90, 91, 30)
    axQ.set_xticks(xticks)
    xtick_labels = [
        f"{abs(val)}°S" if val < 0 else (f"{val}°N" if val > 0 else "적도")
        for val in xticks
    ]
    axQ.set_xticklabels(xtick_labels, fontsize=8)
    axQ.tick_params(axis="y", labelsize=8)
    figQ.tight_layout(pad=0.5)
    st.pyplot(figQ)

with top_col_sky:
    st.subheader("🌤️ 태양의 하늘 경로")
    fig_sky = draw_sun_path(phi_deg, delta, epsilon_deg)
    st.pyplot(fig_sky)

st.subheader("🌞 현재 태양 위치 정보")
st.markdown(info_table_html, unsafe_allow_html=True)

if st.session_state.animate:
    st.session_state.N = (st.session_state.N + 1) % 365
    time.sleep(st.session_state.anim_speed / 1000.0)
    trigger_rerun()
