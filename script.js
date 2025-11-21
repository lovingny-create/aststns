const G = 6.67430e-11;
const SIGMA = 5.670374419e-8;
const SOLAR_RADIUS = 6.957e8; // m
const SOLAR_MASS = 1.98847e30; // kg

const elements = {
  r1: document.getElementById('r1'),
  r2: document.getElementById('r2'),
  t1: document.getElementById('t1'),
  t2: document.getElementById('t2'),
  m1: document.getElementById('m1'),
  m2: document.getElementById('m2'),
  ecc: document.getElementById('ecc'),
  inc: document.getElementById('inc'),
  period: document.getElementById('period'),
};

const labels = {
  r1: document.getElementById('r1-val'),
  r2: document.getElementById('r2-val'),
  t1: document.getElementById('t1-val'),
  t2: document.getElementById('t2-val'),
  m1: document.getElementById('m1-val'),
  m2: document.getElementById('m2-val'),
  ecc: document.getElementById('ecc-val'),
  inc: document.getElementById('inc-val'),
  period: document.getElementById('period-val'),
};

function updateLabels() {
  labels.r1.textContent = Number(elements.r1.value).toFixed(2);
  labels.r2.textContent = Number(elements.r2.value).toFixed(2);
  labels.t1.textContent = `${Number(elements.t1.value).toFixed(0)}`;
  labels.t2.textContent = `${Number(elements.t2.value).toFixed(0)}`;
  labels.m1.textContent = Number(elements.m1.value).toFixed(2);
  labels.m2.textContent = Number(elements.m2.value).toFixed(2);
  labels.ecc.textContent = Number(elements.ecc.value).toFixed(2);
  labels.inc.textContent = `${Number(elements.inc.value).toFixed(1)}°`;
  labels.period.textContent = Number(elements.period.value).toFixed(2);
}

function keplerEccentricAnomaly(meanAnomaly, e) {
  let E = meanAnomaly;
  for (let iter = 0; iter < 12; iter += 1) {
    const f = E - e * Math.sin(E) - meanAnomaly;
    const fPrime = 1 - e * Math.cos(E);
    E -= f / fPrime;
  }
  return E;
}

function trueAnomaly(eccentricAnomaly, e) {
  const sinE = Math.sin(eccentricAnomaly / 2);
  const cosE = Math.cos(eccentricAnomaly / 2);
  return 2 * Math.atan2(Math.sqrt(1 + e) * sinE, Math.sqrt(1 - e) * cosE);
}

function projectedSeparation(a, e, nu, inclination) {
  const r = (a * (1 - e * e)) / (1 + e * Math.cos(nu));
  const x = r * Math.cos(nu);
  const y = r * Math.sin(nu) * Math.cos(inclination);
  const z = r * Math.sin(nu) * Math.sin(inclination);
  return { separation: Math.hypot(x, y), lineOfSight: z };
}

function overlapArea(r1, r2, d) {
  if (d >= r1 + r2) return 0;
  if (d <= Math.abs(r1 - r2)) return Math.PI * Math.min(r1, r2) ** 2;
  const r1Sq = r1 * r1;
  const r2Sq = r2 * r2;
  const alpha = Math.acos((d * d + r1Sq - r2Sq) / (2 * d * r1));
  const beta = Math.acos((d * d + r2Sq - r1Sq) / (2 * d * r2));
  return r1Sq * alpha + r2Sq * beta - d * r1 * Math.sin(alpha);
}

function baselineLuminosity(radius, temperature) {
  const area = Math.PI * radius * radius;
  return SIGMA * area * temperature ** 4;
}

function computeCurve(params) {
  const {
    radius1, radius2, temp1, temp2, mass1, mass2, eccentricity, inclination, periodDays,
  } = params;
  const periodSeconds = periodDays * 86400;
  const a = Math.cbrt(G * (mass1 + mass2) * (periodSeconds / (2 * Math.PI)) ** 2);

  const base1 = baselineLuminosity(radius1, temp1);
  const base2 = baselineLuminosity(radius2, temp2);
  const baseline = base1 + base2;

  const phases = [];
  const normalizedFlux = [];

  const samples = 400;
  for (let idx = 0; idx <= samples; idx += 1) {
    const phase = idx / samples;
    const M = 2 * Math.PI * phase;
    const E = keplerEccentricAnomaly(M, eccentricity);
    const nu = trueAnomaly(E, eccentricity);
    const { separation, lineOfSight } = projectedSeparation(a, eccentricity, nu, inclination);

    const d = separation;
    let flux = baseline;

    const areaOverlap = overlapArea(radius1, radius2, d);
    if (areaOverlap > 0) {
      const star2InFront = lineOfSight >= 0;
      if (star2InFront) {
        const blockedFraction = areaOverlap / (Math.PI * radius1 * radius1);
        flux = base1 * (1 - blockedFraction) + base2;
      } else {
        const blockedFraction = areaOverlap / (Math.PI * radius2 * radius2);
        flux = base1 + base2 * (1 - blockedFraction);
      }
    }

    phases.push(phase);
    normalizedFlux.push(flux / baseline);
  }

  return { phases, normalizedFlux };
}

function redraw() {
  updateLabels();
  const params = {
    radius1: Number(elements.r1.value) * SOLAR_RADIUS,
    radius2: Number(elements.r2.value) * SOLAR_RADIUS,
    temp1: Number(elements.t1.value),
    temp2: Number(elements.t2.value),
    mass1: Number(elements.m1.value) * SOLAR_MASS,
    mass2: Number(elements.m2.value) * SOLAR_MASS,
    eccentricity: Number(elements.ecc.value),
    inclination: (Number(elements.inc.value) * Math.PI) / 180,
    periodDays: Number(elements.period.value),
  };

  const { phases, normalizedFlux } = computeCurve(params);

  const trace = {
    x: phases,
    y: normalizedFlux,
    mode: 'lines',
    line: { color: '#38bdf8', width: 2 },
    hovertemplate: '위상: %{x:.3f}<br>정규화 광도: %{y:.4f}<extra></extra>',
  };

  const layout = {
    margin: { t: 16, r: 12, b: 50, l: 60 },
    xaxis: { title: 'Orbital Phase', range: [0, 1], zeroline: false },
    yaxis: { title: 'Normalized Flux', zeroline: false },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#e5e7eb' },
  };

  Plotly.react('curve', [trace], layout, { responsive: true });
}

Object.values(elements).forEach((el) => {
  el.addEventListener('input', redraw);
});

updateLabels();
redraw();

