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
def mean_anomaly(N: int, M0_deg: float = -3.0) -> float:
    deg_per_day = 360.0 / 365.25
    M_deg = M0_deg + deg_per_day * N
    return math.radians(M_deg)


def eccentric_anomaly(M: float, e: float, n_iter: int = 6) -> float:
    E = M
    for _ in range(n_iter):
        f = E - e * math.sin(E) - M
        fprime = 1 - e * math.cos(E)
        E = E - f / fprime
    return E


def true_anomaly(E: float, e: float) -> float:
    num = math.sqrt(1 + e) * math.sin(E / 2)
    den = math.sqrt(1 - e) * math.cos(E / 2)
    return 2 * math.atan2(num, den)


def solar_declination(v: float, omega_deg: float, epsilon_deg: float) -> Tuple[float, float]:
    lam = v + math.radians(omega_deg) + math.pi  # Earth-Sun phase shift
    eps = math.radians(epsilon_deg)
    delta = math.asin(math.sin(eps) * math.sin(lam))
    return delta, lam


# ============================================
# 2. Insolation
# ============================================
S0 = 1361


def daily_insolation(phi_rad: np.ndarray, delta: float, e: float, lam: float, omega_deg: float) -> np.ndarray:
    rfac = ((1 + e * math.cos(lam - math.radians(omega_deg))) / (1 - e * e)) ** 2

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


# ============================================
# 3. Orbit visualization
# ============================================
def draw_orbit(e: float, omega_deg: float, E_now: float, epsilon_deg: float):

    # 시각적으로 보기 좋게 이심률 강조
    e_vis = min(e * 10, 0.9)

    a = 1.0
    b = a * math.sqrt(1 - e_vis * e_vis)
    omega = math.radians(omega_deg)
    eps = math.radians(epsilon_deg)

    E_all = np.linspace(0, 2 * np.pi, 500)
    x = -a * (np.cos(E_all) - e_vis)
    y = b * np.sin(E_all)

    xR = x * np.cos(omega) - y * np.sin(omega)
    yR = x * np.sin(omega) + y * np.cos(omega)

    xE = -a * (math.cos(E_now) - e_vis)
    yE = b * math.sin(E_now)
    xE_R = xE * np.cos(omega) - yE * np.sin(omega)
    yE_R = xE * np.sin(omega) + yE * np.cos(omega)

    peri_x, peri_y = -a * (1 - e_vis), 0
    ap_x, ap_y = a * (1 + e_vis), 0
    peri_xR = peri_x * np.cos(omega) - peri_y * np.sin(omega)
    peri_yR = peri_x * np.sin(omega) + peri_y * np.cos(omega)
    ap_xR = ap_x * np.cos(omega) - ap_y * np.sin(omega)
    ap_yR = ap_x * np.sin(omega) + ap_y * np.cos(omega)

    fig, ax = plt.subplots(figsize=(5, 5), facecolor="#0a0f1c")
    ax.set_facecolor("#0a0f1c")

    ax.plot(xR, yR, linestyle="--", color="#4db7ff", linewidth=2)
    ax.scatter(0, 0, s=320, color="#0d5c84", edgecolors="white", linewidths=2, label="태양")
    ax.scatter(xE_R, yE_R, s=120, color="#78ffba", edgecolors="#0a0f1c", linewidths=1.5, label="지구")

    # 자전축
    L = 0.35
    dx = L * math.sin(eps)
    dy = L * math.cos(eps)
    ax.plot(
        [xE_R - dx / 2, xE_R + dx / 2],
        [yE_R - dy / 2, yE_R + dy / 2],
        color="white",
        linewidth=2,
    )

    # 근일점/원일점 라벨
    ax.scatter(peri_xR, peri_yR, s=140, color="#78ffba", alpha=0.7)
    ax.scatter(ap_xR, ap_yR, s=140, color="#78ffba", alpha=0.7)
    ax.text(peri_xR - 0.05, peri_yR + 0.18, "근일점", color="white", ha="right", fontsize=11, weight="bold")
    ax.text(ap_xR + 0.05, ap_yR + 0.18, "원일점", color="white", ha="left", fontsize=11, weight="bold")

    ax.set_aspect("equal")
    R = 1 + e_vis + 0.55
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_title("지구 공전 궤도", color="white", fontsize=14)
    ax.tick_params(colors="white", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#1c2f46")
    ax.grid(color="#1c2f46", linestyle="--", linewidth=0.7)
    fig.tight_layout(pad=1.0)

    return fig


# ============================================
# STREAMLIT APP (학생용)
# ============================================
st.set_page_config(layout="wide")
ensure_korean_font()

st.title("밀란코비치 주기에 따른 기후 변화")

if "animate" not in st.session_state:
    st.session_state.animate = False
if "N" not in st.session_state:
    st.session_state.N = 80
if "month" not in st.session_state:
    st.session_state.month = 3
if "day" not in st.session_state:
    st.session_state.day = 21
if "e" not in st.session_state:
    st.session_state.e = 0.0167
if "omega_deg" not in st.session_state:
    st.session_state.omega_deg = 102.9372
if "epsilon_deg" not in st.session_state:
    st.session_state.epsilon_deg = 23.44
if "phi_deg" not in st.session_state:
    st.session_state.phi_deg = 37.0

# --------------------------------------------
# 입력 UI
# --------------------------------------------
with st.sidebar:
    st.header("입력값")

    st.subheader("날짜 선택")

    # 애니메이션 중에는 month/day 자동 갱신
    if st.session_state.animate:
        m, d = month_day_from_day_of_year(st.session_state.N)
        st.session_state.month = m
        st.session_state.day = d

    c_month, c_day = st.columns(2)
    month = int(
        c_month.selectbox("월", list(range(1, 13)), index=st.session_state.month - 1)
    )
    max_day = DAYS_IN_MONTH[month - 1]
    day = int(
        c_day.selectbox(
            "일",
            list(range(1, max_day + 1)),
            index=min(st.session_state.day, max_day) - 1,
        )
    )

    st.session_state.month = month
    st.session_state.day = day

    # 📌 절기 바로가기 버튼
    st.subheader("📌 절기 바로가기")

    cA, cB, cC, cD = st.columns(4)

    if cA.button("춘분"):
        st.session_state.month = 3
        st.session_state.day = 20
        st.session_state.animate = False
        trigger_rerun()

    if cB.button("하지"):
        st.session_state.month = 6
        st.session_state.day = 21
        st.session_state.animate = False
        trigger_rerun()

    if cC.button("추분"):
        st.session_state.month = 9
        st.session_state.day = 22
        st.session_state.animate = False
        trigger_rerun()

    if cD.button("동지"):
        st.session_state.month = 12
        st.session_state.day = 21
        st.session_state.animate = False
        trigger_rerun()

    st.subheader("공전 매개변수")
    e = st.slider("이심률 e", 0.0, 0.1, st.session_state.e, 0.0001, key="e")
    omega_deg = st.slider(
        "세차(ω)", 0.0, 360.0, st.session_state.omega_deg, key="omega_deg"
    )
    epsilon_deg = st.slider(
        "축 경사(ε)", 0.0, 40.0, st.session_state.epsilon_deg, key="epsilon_deg"
    )

    st.subheader("관측자 위도")
    phi_deg = st.slider("위도", -90.0, 90.0, st.session_state.phi_deg, key="phi_deg")

    st.divider()
    anim_speed = st.slider("애니메이션 속도(ms)", 1, 200, 30)

# --------------------------------------------
# 메인 패널
# --------------------------------------------
colL, colR = st.columns([1, 1])

with colL:
    # 날짜 → N 변환
    active_N = (
        day_of_year(st.session_state.month, st.session_state.day)
        if not st.session_state.animate
        else st.session_state.N
    )

    M = mean_anomaly(active_N)
    E_val = eccentric_anomaly(M, e)
    v = true_anomaly(E_val, e)
    delta, lam = solar_declination(v, omega_deg, epsilon_deg)

    fig_orbit = draw_orbit(e, omega_deg, E_val, epsilon_deg)
    st.pyplot(fig_orbit)

with colR:
    st.subheader("🌞 현재 태양 위치 정보")

    st.markdown(f"**입력 날짜:** {st.session_state.month}월 {st.session_state.day}일")
    st.markdown(f"**태양 적위:** {math.degrees(delta):.2f}°")
    st.markdown(f"**위도:** {phi_deg}°")

    # 일사량
    st.subheader("📈 위도별 하루 태양 에너지량 (W/m²)")
    phi_list = np.linspace(-90, 90, 181)
    phi_rad_all = np.radians(phi_list)
    Q = daily_insolation(phi_rad_all, delta, e, lam, omega_deg)

    figQ, axQ = plt.subplots(figsize=(5, 3.2))
    axQ.plot(phi_list, Q)
    axQ.axvline(phi_deg, color="red", linestyle="--")
    axQ.grid(alpha=0.3)
    axQ.set_xlabel("위도")
    axQ.set_ylabel("일사량(W/m²)")
    figQ.tight_layout(pad=0.5)
    st.pyplot(figQ)

    # 남중고도
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

    st.subheader("🌅 남중고도 · 일사량 · 낮 길이")
    c_alt, c_q, c_daylen = st.columns(3)
    c_alt.metric("정오 고도", f"{alpha_noon:.2f}°")
    c_q.metric("일사량", f"{Q_at_lat:.0f} W/m²")
    c_daylen.metric("낮 길이", f"{hours}시간 {minutes:02d}분")


# --------------------------------------------
# 애니메이션 컨트롤
# --------------------------------------------
st.subheader("⏯ 날짜 자동 변화")

c1, c2, c3 = st.columns(3)
if c1.button("▶ Start"):
    st.session_state.N = day_of_year(st.session_state.month, st.session_state.day)
    st.session_state.animate = True
if c2.button("⏸ Pause"):
    st.session_state.animate = False
if c3.button("↺ 1월 1일"):
    st.session_state.N = 0
    st.session_state.animate = False

if st.session_state.animate:
    st.session_state.N = (st.session_state.N + 1) % 365
    time.sleep(anim_speed / 1000.0)
    trigger_rerun()
