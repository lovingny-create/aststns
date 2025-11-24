# Eclipsing binary light-curve simulator

이 프로젝트는 질량, 반지름, 표면 온도, 경사각, 이심률과 같은 기본 물리 파라미터를 입력으로 받아 식쌍성(eclipsing binary)의 광도 곡선을 계산하는 간단한 시뮬레이터입니다. 구면별과 등방성 복사(림브 다크닝 미적용)를 가정하고, 궤도 운동은 케플러 방정식을 통해 처리합니다.

## 구성
- `index.html` / `style.css` / `script.js`: Plotly.js 그래프와 슬라이더 UI를 갖춘 실시간 광도곡선 웹 시뮬레이터
- `binary_lightcurve.py`: 별, 식쌍성 궤도, 광도 곡선 계산 로직을 담은 모듈
- `simulate.py`: 명령행에서 시뮬레이션을 실행하고 CSV·PNG 그래프로 광도 곡선을 저장하는 스크립트
- `app.py`: Streamlit 기반으로 슬라이더를 제공하는 실시간 광도곡선 데스크톱 웹앱
- `proper_motion_parallax.html`: 고유운동·연주시차 시뮬레이터 단일 HTML
- `pages/proper_motion_parallax.py`: 위 HTML을 Streamlit 멀티페이지에서 불러오는 래퍼

## 빠른 시작

### 웹에서 실시간 시뮬레이션
1. 저장소 루트에서 `index.html`을 브라우저로 직접 열거나 정적 서버로 제공합니다.
2. 반지름(r1, r2), 표면온도(T1, T2), 질량(M1, M2), 이심률(e), 경사각(i), 궤도주기(P) 슬라이더를 조정합니다.
3. Plotly.js 그래프가 케플러 방정식(진근점 각 포함)과 원 겹침 면적을 이용해 eclipse flux를 실시간으로 업데이트합니다.
4. 고유운동·연주시차 시뮬레이터는 `proper_motion_parallax.html`을 브라우저로 직접 열면 바로 확인할 수 있으며, `python -m http.server 8000` 후 `http://localhost:8000/proper_motion_parallax.html`로 접속해도 됩니다.

### Streamlit 앱으로 실시간 시뮬레이션
1. 의존성 설치: `pip install streamlit plotly` (matplotlib 없이도 실행 가능)
2. 아래 명령으로 앱 실행: `streamlit run app.py`
3. 기본 페이지에서는 식쌍성 광도곡선 파라미터(반지름, 온도, 경사각, 이심률, 질량, 장반경, 샘플 개수)를 조정하면 그래프와 주기가 즉시 갱신됩니다.
4. 상단 멀티페이지 메뉴에서 **고유운동 + 연주시차 시뮬레이터**를 선택하면 기존 단일 HTML을 그대로 포함한 인터랙티브 도구를 사용할 수 있습니다.

### CLI로 CSV 저장
```bash
python simulate.py --mass1 1.989e30 --mass2 1.5e30 \
  --radius1 6.957e8 --radius2 5e8 \
  --temp1 6000 --temp2 5000 \
  --semi-major-axis 2.5e10 --eccentricity 0.1 \
  --inclination $(python - <<'PY'
import math
print(math.radians(87))
PY
) \
  --phases 400 --output light_curve.csv
```

실행 후 `light_curve.csv`에는 위상(0~1)과 정규화된 상대 광도가 저장되고, 동일한 이름의 PNG(`light_curve.png`)로 광도 곡선 플롯이 생성됩니다. `--figure`로 PNG 경로를 바꿀 수 있으며, 표준 출력에는 궤도 주기(일 단위)도 표시됩니다. 시스템에 matplotlib이 없다면 같은 이름의 SVG로 대체 저장합니다.

## 모델 가정 및 한계
- 별은 완전한 구 형태이며 림브 다크닝, 타이달 변형, 반사 효과 등을 고려하지 않습니다.
- 궤도는 점질량 두 개의 케플러 운동으로 모델링하며, 분광선 속도나 시간 지연 효과는 포함하지 않습니다.
- 광도는 각 별의 표면 밝기를 스테판-볼츠만 법칙으로 계산한 후 식 현상 시 기하학적 가려짐에 따른 면적 감소만을 반영합니다.

## 개발/테스트
```bash
python -m py_compile binary_lightcurve.py simulate.py
```
