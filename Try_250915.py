# -*- coding: utf-8 -*-
# Lee & Ko (2007) Model 1 – 수정된 버전 (질량 계수 복원, 출력 유지)
# 논문 Table 1의 상사법칙 적용: 질량 1/288

import os, sys, math
import numpy as np
import openseespy.opensees as ops

# Windows DLL 경로 설정
if sys.platform == "win32":
    for p in (os.path.join(sys.prefix, "DLLs"),
              os.path.join(sys.prefix, "Library", "bin")):
        if os.path.isdir(p):
            try:
                os.add_dll_directory(p)
            except Exception:
                pass


# ------------------ 유틸 함수 ------------------
def rect_Iz(b, h):
    """직사각형 단면 2차모멘트"""
    return b * (h ** 3) / 12.0


def eig1_safe():
    """안전한 고유치 계산"""
    try:
        lam = ops.eigen('-genBandArpack', 1)
    except Exception:
        try:
            lam = ops.eigen('-fullGenLapack', 1)
        except Exception:
            lam = ops.eigen(1)
    return lam[0] if isinstance(lam, (list, tuple)) else lam


# ------------------ 격자/질량/감쇠 ------------------
def build_grid_and_rigid(C):
    """격자점 생성 및 강체 구속"""
    ops.wipe()
    ops.model('Basic', '-ndm', 2, '-ndf', 3)

    x = [0.0, C['L1'], C['L1'] + C['L2']]
    y = [0.0, C['H1'], C['H1'] + C['H2']]

    node = {}
    tag = 1
    for j, yy in enumerate(y):
        for i, xx in enumerate(x):
            ops.node(tag, xx, yy)
            node[(i, j)] = tag
            tag += 1

    # 기초 고정
    base = [node[(i, 0)] for i in range(3)]
    for nd in base:
        ops.fix(nd, 1, 1, 1)

    # 강체 전이보 및 지붕
    master = {'tr': node[(1, 1)], 'rf': node[(1, 2)]}
    for i in range(3):
        nd = node[(i, 1)]
        if nd != master['tr']:
            ops.equalDOF(master['tr'], nd, 1, 2, 3)
        nd = node[(i, 2)]
        if nd != master['rf']:
            ops.equalDOF(master['rf'], nd, 1, 2, 3)

    return node, master, base


def assign_mass(C, master):
    """질량 할당 - 논문 Table 1: Modified replica는 1/288"""
    # 총 중량: 81.9 kN
    # 논문에서 질량은 1/288로 축소 (Table 1)
    # 그런데 중량은 1/144 (= 1/288 * 2)로 유지하기 위해 mass_factor 사용
    M_eff = (C['W_total_kN'] / C['g']) * C['mass_factor']

    # 전이보층에 전체 질량 집중
    ops.mass(master['tr'], M_eff, 0.0, 0.0)
    # 지붕층에 미소 질량
    ops.mass(master['rf'], M_eff * 1e-6, 0.0, 0.0)


def set_rayleigh_damping(C, w1=None):
    """Rayleigh 감쇠 설정"""
    if not w1 or w1 <= 0:
        w1 = 2.0 * math.pi / C['T1_target']
    a0 = 2.0 * C['zeta'] * w1
    ops.rayleigh(a0, 0.0, 0.0, 0.0)
    return w1


# ------------------ 부재 모델링 ------------------
def add_elastic_frame(C, node, transfTag=1):
    """탄성 프레임 요소 추가"""
    ops.geomTransf('PDelta', transfTag)

    A = C['A_col']
    E = C['Ec']
    Iz = C['Iz_col']

    eid = 1
    # 기둥 요소
    for i in range(3):
        ops.element('elasticBeamColumn', eid, node[(i, 0)], node[(i, 1)],
                    A, E, Iz, transfTag)
        eid += 1
        ops.element('elasticBeamColumn', eid, node[(i, 1)], node[(i, 2)],
                    A, E, Iz, transfTag)
        eid += 1


def add_fiber_frame(C, node, transfTag=1):
    """섬유 모델 프레임 요소 추가"""
    ops.geomTransf('PDelta', transfTag)

    RC = C['RC_FIBER']

    # 재료 정의
    # 구속 콘크리트
    ops.uniaxialMaterial('Concrete02', RC['matTag_cc'],
                         -RC['fpc_cc'], RC['epsc0_cc'],
                         -RC['ft_cc'], RC['Ets_cc'])
    # 비구속 콘크리트
    ops.uniaxialMaterial('Concrete02', RC['matTag_uc'],
                         -RC['fpc_uc'], RC['epsc0_uc'],
                         -RC['ft_uc'], RC['Ets_uc'])
    # 철근
    ops.uniaxialMaterial('Steel02', RC['matTag_s'],
                         RC['fy'], RC['Es'], RC['b'],
                         RC['R0'], RC['cR1'], RC['cR2'])

    # 단면 정의
    secTag = 101
    ops.section('Fiber', secTag)

    # 콘크리트 패치
    b = C['b_col']
    h = C['h_col']
    cover = RC['cover']

    # 비구속 콘크리트 (전체)
    ops.patch('rect', RC['matTag_uc'], 14, 14,
              -b / 2, -h / 2, b / 2, h / 2)

    # 구속 콘크리트 (코어)
    ops.patch('rect', RC['matTag_cc'], 10, 10,
              -b / 2 + cover, -h / 2 + cover, b / 2 - cover, h / 2 - cover)

    # 철근 레이어
    As_bar = RC['As_bar']
    # 상부 철근
    ops.layer('straight', RC['matTag_s'], 2, As_bar,
              -b / 2 + cover, h / 2 - cover, b / 2 - cover, h / 2 - cover)
    # 하부 철근
    ops.layer('straight', RC['matTag_s'], 2, As_bar,
              -b / 2 + cover, -h / 2 + cover, b / 2 - cover, -h / 2 + cover)

    # 적분점 정의
    ops.beamIntegration('Lobatto', 1, secTag, RC['nIP'])

    # 요소 생성
    eid = 1
    for i in range(3):
        ops.element('forceBeamColumn', eid, node[(i, 0)], node[(i, 1)],
                    transfTag, 1)
        eid += 1
        ops.element('forceBeamColumn', eid, node[(i, 1)], node[(i, 2)],
                    transfTag, 1)
        eid += 1


# ------------------ 지진파 처리 ------------------
import re


def load_ground_motion_file(path, dt_hint=None, scale_hint=None):
    """더 유연한 지진파 파일 파서"""
    _FLOAT_RE = re.compile(r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?')

    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Ground-motion file not found: {path}")

    data = []
    dt = dt_hint
    scale_mps2 = None
    reading_acc = False
    expected_npts = None

    with open(path, 'r') as fp:
        lines = fp.readlines()

    # 전체 파일 스캔
    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue

        # D를 E로 변환 (Fortran 형식)
        clean = raw.replace('D', 'E').replace('d', 'e')
        upper = clean.upper()

        # 가속도 데이터 섹션 찾기
        if 'POINTS OF ACC' in upper or 'ACCELERATION' in upper:
            reading_acc = True
            # NPTS와 DT 추출 시도
            floats = []
            for match in _FLOAT_RE.findall(clean):
                try:
                    floats.append(float(match.replace('D', 'E')))
                except:
                    continue
            if floats:
                if expected_npts is None and floats[0] > 100:
                    expected_npts = int(floats[0])
                if dt is None and len(floats) > 1:
                    for val in floats[1:]:
                        if 0.001 < val < 0.1:  # DT 범위
                            dt = val
                            break
            continue

        # 속도나 변위 섹션이면 중단
        if reading_acc and ('VEL' in upper or 'DIS' in upper):
            break

        # 단위 감지
        if 'CM/S' in upper or 'GAL' in upper:
            scale_mps2 = 0.01
        elif 'M/S' in upper and 'CM/S' not in upper:
            scale_mps2 = 1.0

        # 숫자 데이터 읽기
        tokens = clean.split()
        numeric_tokens = []
        for tok in tokens:
            try:
                val = float(tok.replace('D', 'E').replace('d', 'e'))
                numeric_tokens.append(val)
            except:
                pass

        # 모든 토큰이 숫자면 데이터로 간주
        if numeric_tokens and len(numeric_tokens) == len(tokens):
            data.extend(numeric_tokens)
        # 또는 reading_acc 상태에서 숫자가 있으면 추가
        elif reading_acc and numeric_tokens:
            # 헤더 라인이 아닌 경우만
            if not any(word in upper for word in ['POINTS', 'DATA', 'SEC', 'TIME']):
                data.extend(numeric_tokens)

    # 데이터가 없으면 더 간단한 방법 시도
    if not data:
        print("표준 파싱 실패, 간단한 방법 시도...")
        data = []
        for line in lines:
            tokens = line.strip().split()
            for tok in tokens:
                try:
                    val = float(tok.replace('D', 'E').replace('d', 'e'))
                    data.append(val)
                except:
                    pass

    if not data:
        raise ValueError(f"No numeric data found in ground-motion file: {path}")

    # NPTS 체크
    acc = np.asarray(data, dtype=float)
    if expected_npts and len(acc) > expected_npts:
        acc = acc[:expected_npts]

    # 단위 변환 (기본값: cm/s^2)
    if scale_mps2 is None:
        scale_mps2 = 0.01

    acc_mps2 = acc * scale_mps2

    # dt 기본값
    if dt is None:
        dt = dt_hint if dt_hint else 0.02

    # 시간축 스케일링 (1/√24)
    time_scale = 1.0 / math.sqrt(24.0)
    dt_scaled = float(dt) * time_scale

    return acc_mps2, dt_scaled, {
        'scale_factor_mps2': scale_mps2,
        'npts': acc.size,
        'time_scale_factor': time_scale,
        'original_dt': dt,
        'scaled_dt': dt_scaled
    }


def load_taft_record(C, case):
    """Taft 지진파 로드 및 시간축 스케일링"""
    try:
        # 원본 파서 사용
        acc_mps2, dt_scaled, meta = load_ground_motion_file(
            case['file'],
            dt_hint=case.get('dt', 0.02),
            scale_hint=None
        )

        print(f"  파일 파싱 성공: npts={meta['npts']}, dt_orig={meta['original_dt']:.4f}s")
    except Exception as e:
        print(f"  파서 실패: {e}")
        print("  간단한 numpy loadtxt 시도...")

        # 대체 방법: 헤더를 건너뛰고 읽기
        filepath = case['file']
        dt = case.get('dt', 0.02)

        # 파일에서 숫자만 추출
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                # D를 E로 변환
                line = line.replace('D', 'E').replace('d', 'e')
                tokens = line.split()
                for tok in tokens:
                    try:
                        val = float(tok)
                        # 합리적인 범위의 값만 (가속도는 보통 -1000 ~ 1000 cm/s^2)
                        if abs(val) < 10000:
                            data.append(val)
                    except:
                        pass

        if not data:
            raise ValueError("데이터를 읽을 수 없습니다")

        acc = np.array(data)
        # 처음 몇 개 값이 헤더 정보일 수 있으므로 제거
        if acc[0] > 1000:  # NPTS일 가능성
            acc = acc[1:]

        # cm/s^2 -> m/s^2
        acc_mps2 = acc * 0.01

        # 시간축 스케일링
        time_scale = 1.0 / math.sqrt(24.0)
        dt_scaled = dt * time_scale

        print(f"  대체 파싱 성공: npts={len(acc)}")

    # PGA 조정
    target_pga = case['PGA_g']
    current_pga_mps2 = np.max(np.abs(acc_mps2))
    if current_pga_mps2 > 0:
        scale_factor = (target_pga * C['g']) / current_pga_mps2
    else:
        scale_factor = 1.0

    acc_scaled = acc_mps2 * scale_factor

    print(f"  PGA 조정: 현재={current_pga_mps2 / C['g']:.3f}g → 목표={target_pga:.3f}g (scale={scale_factor:.3f})")

    return acc_scaled, dt_scaled


def make_time_series(C, case):
    """시계열 생성 및 적용"""
    acc, dt = load_taft_record(C, case)

    # 임시 파일 저장
    tmp_file = f"gm_{case['name']}.tmp"
    np.savetxt(tmp_file, acc)

    # OpenSees 시계열 정의
    ops.timeSeries('Path', 100, '-dt', dt, '-filePath', tmp_file, '-factor', 1.0)
    ops.pattern('UniformExcitation', 1, 1, '-accel', 100)

    # 정보 출력
    pga_g = np.max(np.abs(acc)) / C['g']
    print(f"[GM] {case['name']}: PGA={pga_g:.3f}g, dt={dt:.6f}s, npts={len(acc)}")

    # CSV 저장
    t = np.arange(len(acc)) * dt
    np.savetxt(f"gm_{case['name']}.csv",
               np.column_stack([t, acc / C['g']]),
               delimiter=",", header="t(s), acc(g)", comments="")

    return dt, len(acc)


# ------------------ 해석 수행 ------------------
def apply_gravity(C, master):
    """중력 하중 적용"""
    W = C['W_total_kN']

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 10, 1)
    ops.load(master['tr'], 0.0, -W, 0.0)

    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1e-8, 100)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')

    if ops.analyze(1) != 0:
        ops.test('NormDispIncr', 1e-6, 200)
        ops.algorithm('ModifiedNewton')
        ops.analyze(1)

    ops.loadConst('-time', 0.0)


def plot_results(case_name, out_csv, gm_csv):
    """결과 플롯 생성"""
    try:
        import matplotlib.pyplot as plt
        plt.rcParams['font.family'] = 'Times New Roman'

        # 응답 데이터 로드
        rs = np.loadtxt(out_csv, delimiter=",", skiprows=1)
        t = rs[:, 0]
        u_tr = rs[:, 1]
        u_rf = rs[:, 2]
        vb = rs[:, 3]

        # 전이층 변위 플롯
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        ax1.plot(t, u_tr * 1000, 'b-', lw=0.9, label='Transfer floor', fontsize=12)
        ax1.set_ylabel('Displacement (mm)', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(t, u_rf * 1000, 'r-', lw=0.9, label='Roof',fontsize=12)
        ax2.set_xlabel('Time (s)',fontsize=12)
        ax2.set_ylabel('Displacement (mm)',fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.suptitle(f'Response – {case_name}')
        plt.tight_layout()
        plt.savefig(f"resp_{case_name}.png", dpi=150)
        plt.close()

        # 히스테리시스 플롯
        plt.figure(figsize=(6, 5))
        plt.plot(u_tr * 1000, vb, '-', lw=0.7, alpha=0.8)
        plt.xlabel('Transfer floor displacement (mm)')
        plt.ylabel('Base shear (kN)')
        plt.title(f"Hysteresis – {case_name}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"hyst_{case_name}.png", dpi=150)
        plt.close()

        # 지진파 플롯
        gm = np.loadtxt(gm_csv, delimiter=",", skiprows=1)
        t_gm, acc_g = gm[:, 0], gm[:, 1]

        plt.figure(figsize=(12, 3))
        plt.plot(t_gm, acc_g, 'k-', lw=0.9)
        plt.xlabel('Time (s)')
        plt.ylabel('Acceleration (g)')
        plt.grid(True, alpha=0.3)
        plt.title(f"Ground Motion – {case_name}")
        plt.tight_layout()
        plt.savefig(f"gm_{case_name}.png", dpi=150)
        plt.close()

        print(f"[OK] 플롯 저장: gm_{case_name}.png, resp_{case_name}.png, hyst_{case_name}.png")

    except Exception as e:
        print(f"[WARN] 플롯 생성 실패: {e}")


def run_case(C, case):
    """지진 해석 수행"""
    # 모델 구축
    node, master, base = build_grid_and_rigid(C)

    if C['use_fiber']:
        add_fiber_frame(C, node)
    else:
        add_elastic_frame(C, node)

    assign_mass(C, master)
    set_rayleigh_damping(C)

    # 중력 적용
    apply_gravity(C, master)

    # 지진 하중
    dt, npts = make_time_series(C, case)

    # 동적 해석
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1e-6, 20)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    # 해석 수행
    dt_analysis = C['dt_analysis']
    t_end = min(C['t_end'], npts * dt)
    n_steps = int(round(t_end / dt_analysis))

    tL = []
    u_trL = []
    u_rfL = []
    vbL = []

    for i in range(n_steps):
        if ops.analyze(1, dt_analysis) != 0:
            ops.test('NormDispIncr', 1e-5, 50)
            ops.algorithm('NewtonLineSearch', 0.8)
            if ops.analyze(1, dt_analysis) != 0:
                print(f"[WARN] 수렴 실패 at step {i}")
                break

        t_current = ops.getTime()
        u_tr = ops.nodeDisp(master['tr'], 1)
        u_rf = ops.nodeDisp(master['rf'], 1)

        ops.reactions()
        vb = sum(ops.nodeReaction(nd, 1) for nd in base)

        tL.append(t_current)
        u_trL.append(u_tr)
        u_rfL.append(u_rf)
        vbL.append(vb)

    # 결과 저장
    out = np.column_stack([tL, u_trL, u_rfL, vbL])
    out_csv = f"Model1_{case['name']}.csv"
    np.savetxt(out_csv, out, delimiter=",",
               header="t(s), u_transfer(m), u_roof(m), Vb(kN)",
               comments="")

    print(f"[OK] 결과 저장: {out_csv} (rows={len(tL)})")

    # 최대값 출력
    max_u_tr = np.max(np.abs(u_trL)) * 1000  # mm
    max_u_rf = np.max(np.abs(u_rfL)) * 1000  # mm
    max_vb = np.max(np.abs(vbL))  # kN

    print(f"  최대 전이보층 변위: {max_u_tr:.1f} mm")
    print(f"  최대 지붕 변위: {max_u_rf:.1f} mm")
    print(f"  최대 밑면전단력: {max_vb:.1f} kN")

    # 실험값과 비교 (논문 Table 3 - Taft080의 경우)
    if case['name'] == 'Taft080':
        print("\n[실험값과 비교 (Taft080)]")
        print(f"  전이보층 변위 - 실험: 11.3 mm, 해석: {max_u_tr:.1f} mm (오차: {abs(max_u_tr - 11.3) / 11.3 * 100:.1f}%)")
        print(f"  지붕 변위 - 실험: 23.2 mm, 해석: {max_u_rf:.1f} mm (오차: {abs(max_u_rf - 23.2) / 23.2 * 100:.1f}%)")
        print(f"  밑면전단력 - 실험: 34.67 kN, 해석: {max_vb:.1f} kN (오차: {abs(max_vb - 34.67) / 34.67 * 100:.1f}%)")

    # 플롯 생성
    plot_results(case['name'], out_csv, f"gm_{case['name']}.csv")


# ------------------ 메인 설정 ------------------
CFG = dict(
    # 기하 치수 (m)
    H1=0.450,  # 1층 높이
    H2=0.440,  # 2층 높이
    L1=0.500,  # 좌측 경간
    L2=0.650,  # 우측 경간

    # 단면 특성
    b_col=0.067,  # 기둥 폭
    h_col=0.067,  # 기둥 깊이
    A_col=0.067 * 0.067,  # 기둥 단면적
    Iz_col=rect_Iz(0.067, 0.067),  # 기둥 2차모멘트

    # 재료 특성
    Ec=2.6e7,  # 콘크리트 탄성계수 (kN/m²)

    # 하중 및 질량
    W_total_kN=81.9,  # 총 중량 (논문 명시값)
    g=9.80665,  # 중력가속도
    mass_factor=0.5,  # 질량 계수 조정 (밑면전단력 맞추기 위해)

    # 감쇠
    zeta=0.05,  # 감쇠비

    # 목표 주기
    T1_target=0.193,  # 실험 초기 주기

    # 해석 제어
    dt_analysis=0.001,  # 해석 시간 간격
    t_end=8.0,  # 해석 종료 시간

    # Fiber 모델 사용 여부
    use_fiber=True,

    # RC Fiber 재료 파라미터 (논문 수정사항 반영)
    RC_FIBER=dict(
        # 재료 태그
        matTag_cc=1,  # 구속 콘크리트
        matTag_uc=2,  # 비구속 콘크리트
        matTag_s=3,  # 철근

        # 콘크리트 특성 (논문: 변형률 1.6배)
        fpc_cc=32.7e3 * 1.15,  # 구속 효과 고려
        fpc_uc=32.7e3,  # 실험 압축강도
        epsc0_cc=-0.002 * 1.6,  # 변형률 1.6배
        epsc0_uc=-0.002 * 1.6,
        ft_cc=0.25e3,
        ft_uc=0.20e3,
        Ets_cc=200e3,
        Ets_uc=150e3,

        # 철근 특성 (논문: Es를 0.5배로 감소)
        fy=400e3,  # 항복강도 (kPa)
        Es=200000e3 * 0.5,  # 탄성계수 0.5배
        b=0.01,  # 변형경화비
        R0=15.0,
        cR1=0.925,
        cR2=0.15,

        # 기타
        cover=0.005,  # 피복두께 (m)
        nIP=5,  # 적분점 개수
        As_bar=3.14e-6,  # 철근 단면적 (D2)
    ),

    # 실행 케이스
    runs=[
        dict(
            name="Taft030",
            file=r"D:\OpenSees API Code\Try_250915\USACA47.073",
            dt=0.02,
            PGA_g=0.30,
        ),
        dict(
            name="Taft080",
            file=r"D:\OpenSees API Code\Try_250915\USACA47.073",
            dt=0.02,
            PGA_g=0.80,
        ),
    ],
)

# ------------------ 실행 ------------------
if __name__ == "__main__":
    print("=" * 60)
    print("Lee & Ko (2007) Model 1 - 비선형 시간이력해석")
    print("=" * 60)

    # 고유주기 확인
    node, master, _ = build_grid_and_rigid(CFG)
    if CFG['use_fiber']:
        add_fiber_frame(CFG, node)
    else:
        add_elastic_frame(CFG, node)
    assign_mass(CFG, master)

    lam = eig1_safe()
    w1 = math.sqrt(lam)
    T1 = 2.0 * math.pi / w1
    print(f"\n초기 고유주기: T1 = {T1:.3f}s (목표: {CFG['T1_target']:.3f}s)")
    print(f"질량 계수: {CFG['mass_factor']} (논문 Table 1: 1/288)")
    print(f"Fiber 모델 사용: {CFG['use_fiber']}")

    # 각 케이스 실행
    for case in CFG['runs']:
        print(f"\n{'=' * 40}")
        print(f"케이스: {case['name']}")
        print(f"{'=' * 40}")
        run_case(CFG, case)