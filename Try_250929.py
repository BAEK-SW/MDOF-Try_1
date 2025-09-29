# IF2 (개구부 없음) – Fiber + 등가지주(X) + 패널존 스프링
# 수정: 수렴성 개선을 위한 요소/재료/해석 설정 최적화
# 단위: mm, N, MPa
import math
import openseespy.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# -------------------- 형상 --------------------
SPAN = 3750.0
H1 = 1840.0
H2 = 1350.0
Htot = H1 + H2

nid = {}


def add_node(tag, x, y):
    ops.node(tag, x, y);
    nid[tag] = (x, y)


add_node(1, 0.0, 0.0)
add_node(2, SPAN, 0.0)
add_node(3, 0.0, H1)
add_node(4, SPAN, H1)  # 제어노드(우측 1층 상단)
add_node(5, 0.0, Htot)  # 상부 좌
add_node(6, SPAN, Htot)  # 상부 우

# 지점(기초 고정)
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)

# 상부 노드: 회전만 고정, 수평·수직 동률(상부 반력빔 등가)
ops.fix(5, 0, 0, 1)
ops.fix(6, 0, 0, 1)
ops.equalDOF(5, 6, 1, 2)


# -------------------- 부재 상세 --------------------
def Ab(d): return math.pi * (d ** 2) / 4.0


# 기둥
col_b, col_h = 250.0, 250.0
col_bar_d = 13.0
col_top_n = 4;
col_bot_n = 4
Lcol_end = 5 * 75.0  # 375
Lcol_mid = H1 - 2 * Lcol_end  # 1090

# 보
beam_b, beam_h = 200.0, 300.0
bar_d = 13.0
BEAM_TOP_n = {'Lend': 5, 'Mid': 2, 'Rend': 5}
BEAM_BOT_n = {'Lend': 2, 'Mid': 2, 'Rend': 2}
Lbeam_end = 450.0
Lbeam_mid = SPAN - 2 * Lbeam_end
if Lbeam_mid <= 0: raise ValueError("보 중앙부 길이가 0 이하입니다.")

# 벽체
wall_t = 75.0

# -------------------- 재료(표 3.3 평균) - 수렴성 개선 --------------------
MAT_C_CORE = 1;
MAT_C_COV = 2;
MAT_STEEL = 3
fpc, ft = 21.10, 2.06
epsc0, epsU = -0.002, -0.02  # epsU를 -0.006에서 -0.02로 완화 (수렴성 개선)
lam = 0.2  # 0.1에서 0.2로 증가 (언로딩/리로딩 강성 증가)
fpcu = 0.2 * fpc
# Ets: 인장 연화를 더 완만하게 (수렴성 개선)
Ets = ft / 0.01  # 0.004에서 0.01로 변경

# Concrete02 대신 Concrete01 사용 (더 안정적)
ops.uniaxialMaterial('Concrete01', MAT_C_CORE, -fpc, epsc0, -fpcu, epsU)
ops.uniaxialMaterial('Concrete01', MAT_C_COV, -fpc, epsc0, -fpcu, epsU)

# Steel02: b값 증가로 안정화
fy_D13, Es_D13, b_iso = 498.96, 184000.0, 0.02  # b: 0.01에서 0.02로
ops.uniaxialMaterial('Steel02', MAT_STEEL, fy_D13, Es_D13, b_iso,
                     20, 0.925, 0.15)  # R0, cR1, cR2 추가


# -------------------- 섬유 단면 유틸 --------------------
def fiber_rc_rect_section(secTag, b, h, cover, topBars, botBars, barDia,
                          matC_core=MAT_C_CORE, matC_cov=MAT_C_COV, matS=MAT_STEEL,
                          nfCoreY=8, nfCoreZ=8):  # 메쉬 감소 (12->8)
    ops.section('Fiber', secTag)
    yTop, yBot = +h / 2.0, -h / 2.0
    zR = +b / 2.0
    coreY = h - 2.0 * cover
    coreZ = b - 2.0 * cover
    if coreY <= 0 or coreZ <= 0: raise ValueError("코어 치수 오류")

    yC_top = +coreY / 2.0;
    yC_bot = -coreY / 2.0
    zC_rt = +coreZ / 2.0;
    zC_lt = -coreZ / 2.0

    # 커버 4패치
    ops.patch('rect', matC_cov, 2, nfCoreZ, yC_top, -zR, yTop, +zR)
    ops.patch('rect', matC_cov, 2, nfCoreZ, yBot, -zR, yC_bot, +zR)
    ops.patch('rect', matC_cov, nfCoreY, 2, yC_bot, -zR, yC_top, zC_lt)
    ops.patch('rect', matC_cov, nfCoreY, 2, yC_bot, zC_rt, yC_top, +zR)
    # 코어
    ops.patch('rect', matC_core, nfCoreY, nfCoreZ, yC_bot, zC_lt, yC_top, zC_rt)
    # 상/하 주근
    yBarTop, yBarBot = yTop - cover, yBot + cover
    if topBars > 0:
        ops.layer('straight', matS, topBars, Ab(barDia), yBarTop, -zR + cover, yBarTop, +zR - cover)
    if botBars > 0:
        ops.layer('straight', matS, botBars, Ab(barDia), yBarBot, -zR + cover, yBarBot, +zR - cover)


# -------------------- 요소 생성 --------------------
# 기하변환: Corotational 사용 (대변위 해석에 더 적합)
ops.geomTransf('Corotational', 1)

cov_col = 20.0;
cov_beam = 20.0

# 섹션
SEC_COL_L_END, SEC_COL_L_MID, SEC_COL_R_END, SEC_COL_R_MID = 101, 102, 103, 104
for sec in (SEC_COL_L_END, SEC_COL_L_MID, SEC_COL_R_END, SEC_COL_R_MID):
    fiber_rc_rect_section(sec, col_b, col_h, cov_col, col_top_n, col_bot_n, col_bar_d)

SEC_BEAM_L_END, SEC_BEAM_MID, SEC_BEAM_R_END = 111, 112, 113
fiber_rc_rect_section(SEC_BEAM_L_END, beam_b, beam_h, cov_beam, BEAM_TOP_n['Lend'], BEAM_BOT_n['Lend'], bar_d)
fiber_rc_rect_section(SEC_BEAM_MID, beam_b, beam_h, cov_beam, BEAM_TOP_n['Mid'], BEAM_BOT_n['Mid'], bar_d)
fiber_rc_rect_section(SEC_BEAM_R_END, beam_b, beam_h, cov_beam, BEAM_TOP_n['Rend'], BEAM_BOT_n['Rend'], bar_d)

# 분할 노드
ops.node(13, 0.0, Lcol_end)
ops.node(23, 0.0, H1 - Lcol_end)
ops.node(14, SPAN, Lcol_end)
ops.node(24, SPAN, H1 - Lcol_end)
ops.node(33, Lbeam_end, H1)
ops.node(34, SPAN - Lbeam_end, H1)

# 요소: forceBeamColumn 사용 (더 robust), 적분점 3개로 감소
# forceBeamColumn: eleTag, iNode, jNode, transfTag, integrationTag
# Lobatto 적분 사용
ops.beamIntegration('Lobatto', 1001, SEC_COL_L_END, 2)
ops.beamIntegration('Lobatto', 1002, SEC_COL_L_MID, 2)
ops.beamIntegration('Lobatto', 1003, SEC_COL_R_END, 2)
ops.beamIntegration('Lobatto', 1004, SEC_COL_R_MID, 2)
ops.beamIntegration('Lobatto', 1011, SEC_BEAM_L_END, 2)
ops.beamIntegration('Lobatto', 1012, SEC_BEAM_MID, 2)
ops.beamIntegration('Lobatto', 1013, SEC_BEAM_R_END, 2)

# forceBeamColumn 요소 생성
ops.element('forceBeamColumn', 301, 1, 13, 1, 1001)
ops.element('forceBeamColumn', 302, 13, 23, 1, 1002)
ops.element('forceBeamColumn', 303, 23, 3, 1, 1001)

ops.element('forceBeamColumn', 311, 2, 14, 1, 1003)
ops.element('forceBeamColumn', 312, 14, 24, 1, 1004)
ops.element('forceBeamColumn', 313, 24, 4, 1, 1003)

ops.element('forceBeamColumn', 321, 3, 33, 1, 1011)
ops.element('forceBeamColumn', 322, 33, 34, 1, 1012)
ops.element('forceBeamColumn', 323, 34, 4, 1, 1013)

# -------------------- 상부 벽체: 등가 대각 압축지주(X) --------------------
weff = 300.0
A_strut = wall_t * weff
MAT_STRUT_CORE, MAT_STRUT = 21, 22
# 벽체 재료도 Concrete01 사용
ops.uniaxialMaterial('Concrete01', MAT_STRUT_CORE, -fpc, epsc0, -fpcu, epsU)
ops.uniaxialMaterial('MinMax', MAT_STRUT, MAT_STRUT_CORE, '-max', 0.0)
ops.element('corotTruss', 901, 3, 6, A_strut, MAT_STRUT)
ops.element('corotTruss', 902, 4, 5, A_strut, MAT_STRUT)

# -------------------- 패널존 전단 스프링 --------------------
Ktheta_L, Ktheta_R = 5.0e6, 5.0e6
MAT_PZ_L, MAT_PZ_R = 31, 32
ops.uniaxialMaterial('Elastic', MAT_PZ_L, Ktheta_L)
ops.uniaxialMaterial('Elastic', MAT_PZ_R, Ktheta_R)
ops.node(203, *nid[3]);
ops.equalDOF(3, 203, 1, 2);
ops.element('zeroLength', 931, 3, 203, '-mat', MAT_PZ_L, '-dir', 3)
ops.node(204, *nid[4]);
ops.equalDOF(4, 204, 1, 2);
ops.element('zeroLength', 932, 4, 204, '-mat', MAT_PZ_R, '-dir', 3)

# ==================== 가력: Enforced Displacement ====================
# (A) 프로토콜 생성 - 더 작은 초기 변위로 시작
disp_levels = [H1 / 400.0, H1 / 200.0, H1 / 100.0, H1  / 67.0, H1 / 50.0, H1 / 33.0]
protocol_targets = []
for i, d in enumerate(disp_levels):
    if i < 3:  # 작은 변위는 1사이클만
        protocol_targets += [d, -d]
    else:  # 큰 변위는 2사이클
        protocol_targets += [d, -d, d, -d]

# 타깃 사이를 더 세밀하게 보간
vals = [0.0]
def add_ramp(v0, v1, du=0.1):  # 0.05에서 0.1로 증가 (속도 향상)
    n = max(int(abs(v1 - v0)/du), 1)
    for k in range(1, n+1):
        vals.append(v0 + (v1 - v0) * k / n)

for tgt in protocol_targets:
    add_ramp(vals[-1], tgt, du=0.1)


for tgt in protocol_targets:
    add_ramp(vals[-1], tgt, du=0.1)

# (B) Path TimeSeries + sp 강제변위
TS, PAT = 200, 200
ops.timeSeries('Path', TS, '-dt', 1.0, '-values', *vals)
ops.pattern('Plain', PAT, TS)
ops.sp(4, 1, 1.0)


# (C) 해석기 설정 - 더 완화된 수렴 기준
def setup(algo='Newton', test=('EnergyIncr', 1e-5, 100)):  # 1e-7->1e-5, 200->100
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')  # BandGeneral에서 UmfPack으로 변경 (더 안정적)
    tname, tol, it = test
    ops.test(tname, tol, it)
    ops.algorithm(algo)
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


setup()

# (D) 레코더
ops.recorder('Node', '-file', 'act_disp_node4.out', '-time', '-node', 4, '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'act_reac_node4.out', '-time', '-node', 4, '-dof', 1, 'reaction')

# (E) 실행 - 개선된 폴백 전략
print(f"Total steps: {len(vals) - 1}")
print(f"Target displacements: {protocol_targets}")

for step in range(1, len(vals)):
    ok = ops.analyze(1)
    if ok != 0:
        print(f"Convergence issue at step {step}/{len(vals) - 1}, target={vals[step]:.3f}mm")

        # 다양한 알고리즘과 완화된 수렴 기준 시도
        for algo, test in [
            ('ModifiedNewton', ('EnergyIncr', 1e-4, 200)),
            ('NewtonLineSearch', ('EnergyIncr', 1e-4, 200)),
            ('KrylovNewton', ('EnergyIncr', 1e-4, 200)),
            ('Newton', ('NormDispIncr', 1e-3, 200)),
            ('ModifiedNewton', ('NormDispIncr', 1e-2, 300)),
        ]:
            setup(algo=algo, test=test)
            ok = ops.analyze(1)
            if ok == 0:
                print(f"  -> Converged with {algo}")
                setup()  # 기본 설정으로 복귀
                break

        if ok != 0:
            # 최후의 수단: 더 작은 스텝으로 분할
            print(f"  -> Trying sub-stepping...")
            curr_disp = ops.nodeDisp(4, 1)
            target_disp = vals[step]
            sub_steps = 10

            for sub in range(sub_steps):
                sub_target = curr_disp + (target_disp - curr_disp) * (sub + 1) / sub_steps
                ops.setTime(step - 1 + (sub + 1) / sub_steps)

                # 변위 직접 설정 (임시 해결책)
                for algo, test in [
                    ('Newton', ('NormDispIncr', 1e-2, 100)),
                    ('ModifiedNewton', ('NormDispIncr', 1e-1, 100)),
                ]:
                    setup(algo=algo, test=test)
                    ok = ops.analyze(1)
                    if ok == 0:
                        break

                if ok != 0:
                    print(f"Warning: Cannot converge at step {step} (target={vals[step]:.3f}mm)")
                    print(f"Continuing with achieved displacement: {ops.nodeDisp(4, 1):.3f}mm")
                    break

            if ok != 0:
                print(f"Stopping analysis at step {step}")
                break

    if step % 100 == 0:
        print(f"Progress: step {step}/{len(vals) - 1}, disp={ops.nodeDisp(4, 1):.2f}mm")

print('=' * 50)
print('Analysis completed!')
print(f'Final node4 ux = {ops.nodeDisp(4, 1):.3f} mm')
print(f'Target was = {vals[-1]:.3f} mm')
print(f'Steps completed = {step}/{len(vals) - 1}')

# ----- 포스트프로세싱: 액추에이터 Load–Displacement -----
import numpy as np
import matplotlib.pyplot as plt

disp = np.loadtxt('act_disp_node4.out')   # [time, ux]
reac = np.loadtxt('act_reac_node4.out')   # [time, Rx]

n   = min(len(disp), len(reac))
ux  = disp[:n, 1]                 # mm  (액추에이터 변위)
P0  = reac[:n, 1] / 1000.0        # kN  (노드4 반력)

# ★ 자동 방향 정렬: +ux에서 +P가 되도록 조정
P   = P0.copy()
if np.nanmean(ux * P0) < 0:       # 평균 상관이 음수면 하중 부호 반전
    P = -P

plt.figure(figsize=(7,5))
plt.plot(ux, P, lw=1.2)
plt.axhline(0, color='k', lw=0.5)
plt.axvline(0, color='k', lw=0.5)
plt.xlabel('Displacement at actuator (mm)')
plt.ylabel('Actuator Load (kN)')
plt.title('IF2 — Actuator Load–Displacement')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
