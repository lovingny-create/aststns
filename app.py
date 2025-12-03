"""Streamlit UI for interactive eclipsing binary light-curve exploration."""
from __future__ import annotations

import math

import plotly.graph_objects as go
import streamlit as st

from simulate import compute_light_curve

SOLAR_RADIUS = 6.957e8
SOLAR_MASS = 1.989e30


def main() -> None:
    st.set_page_config(page_title="Eclipsing Binary Light Curve", layout="wide")
    st.title("식쌍성 광도곡선 시뮬레이터")

    st.markdown(
        """
        반지름, 온도, 경사각, 이심률을 조절하면 실시간으로 정규화 광도곡선과 시선속도 곡선이 함께 업데이트됩니다.
        기본 단위는 태양 단위(반지름, 질량)와 켈빈(온도)이며, 경사각은 도 단위 슬라이더로 설정됩니다.
        """
    )

    with st.sidebar:
        st.header("파라미터")
        radius1 = st.slider("주성 반지름 (R☉)", 0.1, 5.0, 1.0, 0.01)
        radius2 = st.slider("반성 반지름 (R☉)", 0.1, 5.0, 0.7, 0.01)
        temp1 = st.slider("주성 표면온도 (K)", 3000, 15000, 6000, 50)
        temp2 = st.slider("반성 표면온도 (K)", 3000, 15000, 5000, 50)
        inclination_deg = st.slider("경사각 (deg)", 60.0, 90.0, 87.0, 0.1)
        eccentricity = st.slider("이심률", 0.0, 0.9, 0.1, 0.01)
        semi_major_axis = st.number_input(
            "궤도 장반경 (m)",
            min_value=1e9,
            max_value=1e11,
            value=2.5e10,
            step=1e9,
            format="%.2e",
        )
        mass1 = st.number_input(
            "주성 질량 (kg)",
            min_value=0.1 * SOLAR_MASS,
            max_value=10 * SOLAR_MASS,
            value=1.0 * SOLAR_MASS,
            step=0.1 * SOLAR_MASS,
            format="%.2e",
        )
        mass2 = st.number_input(
            "반성 질량 (kg)",
            min_value=0.1 * SOLAR_MASS,
            max_value=10 * SOLAR_MASS,
            value=0.75 * SOLAR_MASS,
            step=0.1 * SOLAR_MASS,
            format="%.2e",
        )
        phases = st.slider("샘플 개수", 100, 800, 400, 50)
        x_axis = st.radio("X축", ["위상", "시간 (일)", "시간 (분)"], index=0)

    phases_values, fluxes, rv_primary, rv_secondary, period = compute_light_curve(
        mass1=mass1,
        mass2=mass2,
        radius1=radius1 * SOLAR_RADIUS,
        radius2=radius2 * SOLAR_RADIUS,
        temp1=float(temp1),
        temp2=float(temp2),
        semi_major_axis=float(semi_major_axis),
        eccentricity=float(eccentricity),
        inclination=math.radians(float(inclination_deg)),
        phases=int(phases),
    )

    if x_axis == "시간 (일)":
        x_values = [phase * period / 86400 for phase in phases_values]
        x_label = "시간 (일)"
        hover_label = "일"
    elif x_axis == "시간 (분)":
        x_values = [phase * period / 60 for phase in phases_values]
        x_label = "시간 (분)"
        hover_label = "분"
    else:
        x_values = phases_values
        x_label = "궤도 위상"
        hover_label = "위상"

    fig = go.Figure()
    flux_hover = (
        "{label}: %{{x:.3f}}<br>정규화 광도: %{{y:.4f}}<extra></extra>"
    ).format(label=hover_label)

    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=fluxes,
            mode="lines",
            line=dict(color="royalblue", width=2),
            name="광도",
            hovertemplate=flux_hover,
        )
    )
    fig.update_layout(
        title="정규화 광도곡선",
        xaxis_title=x_label,
        yaxis_title="상대 광도",
        margin=dict(l=40, r=20, t=50, b=40),
        height=380,
    )

    rv_fig = go.Figure()
    rv_primary_hover = (
        "{label}: %{{x:.3f}}<br>v<sub>r1</sub>: %{{y:.2f}} m/s<extra></extra>"
    ).format(label=hover_label)
    rv_secondary_hover = (
        "{label}: %{{x:.3f}}<br>v<sub>r2</sub>: %{{y:.2f}} m/s<extra></extra>"
    ).format(label=hover_label)

    rv_fig.add_trace(
        go.Scatter(
            x=x_values,
            y=rv_primary,
            mode="lines",
            line=dict(color="indianred", width=2),
            name="주성 RV",
            hovertemplate=rv_primary_hover,
        )
    )
    rv_fig.add_trace(
        go.Scatter(
            x=x_values,
            y=rv_secondary,
            mode="lines",
            line=dict(color="seagreen", width=2),
            name="반성 RV",
            hovertemplate=rv_secondary_hover,
        )
    )
    rv_fig.update_layout(
        title="시선속도 곡선",
        xaxis_title=x_label,
        yaxis_title="시선속도 (m/s)",
        margin=dict(l=40, r=20, t=50, b=40),
        height=380,
    )

    st.plotly_chart(fig, use_container_width=True)
    st.plotly_chart(rv_fig, use_container_width=True)
    st.info(f"궤도 주기: {period / 86400:.2f} 일")


if __name__ == "__main__":
    main()
