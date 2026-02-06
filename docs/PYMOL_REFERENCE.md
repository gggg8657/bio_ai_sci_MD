# PyMOL 완전 가이드

> **Molecular Graphics System**  
> 버전: 3.1.0 (Open-Source)  
> 공식 사이트: https://pymol.org  
> 문서: https://pymol.org/dokuwiki/

---

## 목차

1. [개요](#개요)
2. [설치](#설치)
3. [기본 사용법](#기본-사용법)
4. [표현 방식 (Representations)](#표현-방식-representations)
5. [선택 문법 (Selections)](#선택-문법-selections)
6. [색상 및 렌더링](#색상-및-렌더링)
7. [뷰 제어](#뷰-제어)
8. [정렬 및 중첩](#정렬-및-중첩)
9. [측정](#측정)
10. [분자 편집](#분자-편집)
11. [무비 제작](#무비-제작)
12. [스크립팅](#스크립팅)
13. [고급 기능](#고급-기능)
14. [명령어 레퍼런스](#명령어-레퍼런스)

---

## 개요

PyMOL은 단백질, 핵산, 소분자 등 분자 구조를 시각화하는 **업계 표준** 도구입니다.

### 핵심 특징

| 특징 | 설명 |
|------|------|
| **고품질 렌더링** | 레이트레이싱으로 출판 품질 이미지 |
| **다양한 표현** | Cartoon, Surface, Sticks, Spheres 등 |
| **강력한 선택** | 복잡한 원자 선택 문법 |
| **Python 통합** | 완전한 Python API |
| **구조 정렬** | 서열/구조 기반 중첩 |
| **무비 제작** | 애니메이션 및 비디오 생성 |
| **확장성** | 플러그인 및 스크립트 지원 |

---

## 설치

### Conda (권장)

```bash
conda install -c conda-forge pymol-open-source
```

### 확인

```bash
pymol --version
# PyMOL 3.1.0
```

### WSL에서 실행

WSL에서 GUI 사용 시:
- **Windows 11**: WSLg 자동 지원
- **Windows 10**: X 서버 필요 (VcXsrv, X410 등)

```bash
# X 서버 설정 (Windows 10)
export DISPLAY=:0

# 실행
pymol structure.pdb
```

---

## 기본 사용법

### 명령줄 실행

```bash
# GUI 모드
pymol structure.pdb

# 여러 파일
pymol model1.pdb model2.pdb model3.pdb

# 헤드리스 모드 (이미지만 저장)
pymol -c structure.pdb -d "png output.png; quit"
```

### 기본 명령어

```python
# 파일 로드
load structure.pdb

# 저장
save output.pdb, selection
save image.png

# 종료
quit
```

### 도움말

```python
help                    # 전체 도움말
help commands           # 명령어 목록
help <command>          # 특정 명령어 도움말
help selections         # 선택 문법
```

---

## 표현 방식 (Representations)

### 사용 가능한 표현

| 표현 | 설명 | 용도 |
|------|------|------|
| **lines** | 결합을 선으로 표시 | 기본, 빠른 탐색 |
| **sticks** | 두꺼운 막대 | 리간드, 활성 부위 |
| **spheres** | 원자를 구로 표시 | 공간 충전 모델 |
| **cartoon** | 2차 구조 리본 | 단백질 전체 구조 |
| **ribbon** | 간단한 리본 | 빠른 개요 |
| **surface** | 분자 표면 | 결합 포켓, 상호작용 |
| **mesh** | 메쉬 표면 | 반투명 표면 |
| **dots** | 점 표면 | 가벼운 표면 |
| **labels** | 텍스트 라벨 | 잔기/원자 표시 |
| **nonbonded** | 비결합 원자 | 물, 이온 |
| **nb_spheres** | 비결합 구 | 물 분자 |

### 명령어

```python
# 표현 켜기
show cartoon               # 전체에 cartoon
show sticks, resi 100      # 잔기 100에 sticks
show surface, chain A      # 체인 A에 surface

# 표현 끄기
hide everything            # 모든 표현 끄기
hide lines                 # lines만 끄기
hide cartoon, chain B      # 체인 B의 cartoon 끄기

# 표현 교체 (as = hide all + show)
as cartoon                 # cartoon만 표시
as sticks, organic         # 유기 분자를 sticks로

# 모든 표현 보기
show everything
```

### Cartoon 스타일

```python
# Cartoon 타입 설정
cartoon loop               # 루프 스타일
cartoon tube               # 튜브 스타일
cartoon automatic          # 자동 (기본)
cartoon oval               # 타원형 나선
cartoon rectangle          # 직사각형 시트

# 설정
set cartoon_oval_length, 1.2
set cartoon_rect_length, 1.4
set cartoon_loop_radius, 0.3
```

---

## 선택 문법 (Selections)

### 기본 선택자

| 선택자 | 단축형 | 예시 | 설명 |
|--------|--------|------|------|
| `name` | `n.` | `name CA` | 원자 이름 |
| `resn` | `r.` | `resn ALA` | 잔기 이름 |
| `resi` | `i.` | `resi 100` | 잔기 번호 |
| `chain` | `c.` | `chain A` | 체인 ID |
| `segi` | `s.` | `segi SEG1` | 세그먼트 |
| `elem` | `e.` | `elem C` | 원소 |
| `alt` | - | `alt A` | 대체 위치 |
| `id` | - | `id 100` | 원자 ID |

### 특수 선택자

| 선택자 | 설명 |
|--------|------|
| `all` / `*` | 모든 원자 |
| `none` | 선택 없음 |
| `hydrogen` / `h.` | 수소 원자 |
| `hetatm` | HETATM 레코드 |
| `organic` | 유기 분자 (리간드) |
| `solvent` | 용매 (물) |
| `polymer` | 고분자 (단백질, 핵산) |
| `visible` / `v.` | 보이는 원자 |
| `enabled` | 활성화된 객체 |

### 논리 연산자

```python
# AND
select mysel, chain A and resi 100
select mysel, c. A & i. 100

# OR
select mysel, resn ALA or resn GLY
select mysel, r. ALA | r. GLY

# NOT
select mysel, not hydrogen
select mysel, ! h.

# 복합
select active_site, chain A and resi 100-120 and not hydrogen
```

### 거리 기반 선택

```python
# around: 주변 원자
select near_ligand, organic around 5     # 리간드 5Å 이내

# within: 거리 내 원자
select contacts, chain A within 4 of chain B

# expand: 선택 확장
select expanded, resi 100 expand 3

# gap: 갭 거리
select far, chain A gap 10
```

### 구조 기반 선택

```python
# byres: 잔기 단위로 확장
select whole_res, byres (name CA within 5 of organic)
byres organic around 4

# byobj: 객체 단위로 확장
select whole_obj, byobj chain A

# bycalpha: CA 기준
select backbone, bycalpha all
```

### 범위 및 패턴

```python
# 잔기 범위
resi 10-50                # 10부터 50
resi 10+20+30            # 10, 20, 30

# 체인 여러 개
chain A+B+C

# 와일드카드
resn AL*                 # ALA, ALN 등
name C*                  # C, CA, CB 등
```

### 속성 비교

```python
# B-factor
select high_b, b > 50
select low_b, b < 20

# Occupancy
select partial, q < 1.0

# Formal charge
select positive, formal_charge > 0

# Partial charge
select polar, partial_charge > 0.3 or partial_charge < -0.3
```

### 선택 저장 및 사용

```python
# 선택 생성
select active_site, resi 100-120 and chain A

# 선택 사용
show sticks, active_site
color red, active_site

# 선택 삭제
delete active_site
```

---

## 색상 및 렌더링

### 기본 색상

```python
# 단일 색상
color red, chain A
color blue, organic
color yellow, resi 100

# 원소별 색상 (CPK)
util.cbaw               # 탄소=흰색
util.cbag               # 탄소=회색
util.cbac               # 탄소=청록색
util.cbam               # 탄소=마젠타
util.cbay               # 탄소=노랑
util.cbas               # 탄소=연어색
util.cbap               # 탄소=분홍
```

### 스펙트럼 색상

```python
# 잔기 번호로 스펙트럼
spectrum count, rainbow, chain A

# B-factor로 색상
spectrum b, blue_white_red

# 체인별 색상
util.cbc                # 체인별 색상

# 2차 구조별
util.cbss               # helix=red, sheet=yellow, loop=green
```

### 사용자 정의 색상

```python
# 새 색상 정의
set_color mycolor, [0.8, 0.2, 0.5]
color mycolor, selection

# RGB 직접
color 0x00FF00, selection      # 녹색
```

### 배경 색상

```python
bg_color white
bg_color black
bg_color gray
```

### 렌더링

```python
# 레이트레이싱 (고품질)
ray                            # 현재 크기
ray 1920, 1080                 # 지정 크기
ray 4000, 3000                 # 출판용 고해상도

# 빠른 그리기
draw                           # 안티앨리어싱 없음
draw 1920, 1080

# PNG 저장
png output.png                 # 현재 뷰
png output.png, ray=1          # 레이트레이싱 후 저장
png output.png, 1920, 1080, ray=1
```

### 렌더링 설정

```python
# 안티앨리어싱
set antialias, 2

# 그림자
set ray_shadows, on
set ray_shadow_decay_factor, 0.1

# 조명
set light_count, 4
set spec_reflect, 1.5

# 투명도
set transparency, 0.5, surface
set cartoon_transparency, 0.3
```

---

## 뷰 제어

### 기본 뷰 조작

```python
# 줌
zoom                     # 전체 보기
zoom chain A             # 선택 영역에 맞춤
zoom resi 100, 5         # 5Å 버퍼로 줌

# 중심
center resi 100
origin resi 100          # 회전 중심 설정

# 방향 설정
orient                   # 주축 정렬
orient chain A
```

### 회전 및 이동

```python
# 회전 (도 단위)
turn x, 90               # X축으로 90도
turn y, 45
turn z, 180

# 이동
move x, 10               # X 방향 10단위
move y, -5

# 클리핑
clip near, -5            # 앞 클리핑 조절
clip far, 10
clip slab, 20            # 슬래브 두께
```

### 뷰 저장/복원

```python
# 뷰 저장
view myview, store

# 뷰 복원
view myview, recall

# 현재 뷰 행렬
get_view                 # 18개 숫자 출력
set_view ([...])         # 뷰 설정
```

### Scene 관리

```python
# Scene 저장 (뷰 + 표현 + 색상)
scene F1, store          # F1 키에 저장
scene myname, store

# Scene 복원
scene F1, recall
scene myname

# Scene 목록
scene                    # 목록 출력

# Scene 삭제
scene myname, clear
```

### 창 설정

```python
# 뷰포트 크기
viewport 1920, 1080

# 전체 화면
full_screen on
full_screen off
```

---

## 정렬 및 중첩

### align (서열 기반)

```python
# 기본 정렬
align mobile, target

# 결과 객체 생성
align protA, protB, object=alignment

# 파라미터
align protA, protB, cycles=5, cutoff=2.0

# CA 원자만
align protA////CA, protB////CA
```

### super (구조 기반)

```python
# 서열 유사성 낮을 때 권장
super mobile, target

# 파라미터
super protA, protB, cycles=5, cutoff=2.0, object=super_aln
```

### cealign (CE 알고리즘)

```python
# 가장 강력한 구조 정렬
cealign target, mobile

# 결과
cealign protB, protA, object=ce_aln
```

### pair_fit (원자 쌍 피팅)

```python
# 특정 원자 쌍으로 정렬
pair_fit \
  protA///10/CA, protB///10/CA, \
  protA///20/CA, protB///20/CA, \
  protA///30/CA, protB///30/CA
```

### fit (좌표 피팅)

```python
# 동일 원자 수 필요
fit mobile, target
```

### RMSD 계산

```python
# 현재 좌표 RMSD
rms_cur mobile, target

# 피팅 후 RMSD
rms mobile, target

# 내부 상태 RMSD
intra_rms mobile
intra_fit mobile        # 상태간 피팅
```

---

## 측정

### 거리 측정

```python
# 두 원자 간 거리
distance dist1, /protA/A/100/CA, /protA/A/200/CA

# 선택 간 거리
distance polar_contacts, chain A, chain B, mode=2

# 모드
# mode=0: 모든 원자 쌍
# mode=2: 극성 접촉
# mode=4: 수소 결합
```

### 각도 측정

```python
# 세 원자 각도
angle ang1, \
  /prot/A/100/N, \
  /prot/A/100/CA, \
  /prot/A/100/C
```

### 이면각 측정

```python
# 네 원자 이면각
dihedral dih1, \
  /prot/A/100/N, \
  /prot/A/100/CA, \
  /prot/A/100/C, \
  /prot/A/101/N
```

### 자동 측정 표시

```python
# 극성 접촉 표시
distance hbonds, chain A, chain B, mode=2
show dashes, hbonds

# 스타일 설정
set dash_color, yellow, hbonds
set dash_width, 2.0
set dash_gap, 0.3
```

---

## 분자 편집

### 원자/잔기 조작

```python
# 삭제
remove hydrogen          # 수소 제거
remove solvent           # 물 제거
remove resi 100          # 잔기 삭제

# 이름 변경
set_name old_name, new_name

# 복사
create new_obj, selection

# 추출 (원본에서 제거)
extract new_obj, selection
```

### 구조 수정

```python
# 수소 추가
h_add selection
h_fill                   # 불완전한 원자에 H 추가

# 결합 수정
bond atom1, atom2        # 결합 추가
unbond atom1, atom2      # 결합 제거

# 원자 치환
replace O, N, 3          # N으로 3개 결합

# 회전
rotate x, 45, selection
rotate y, 90, selection, origin=[0,0,0]

# 이동
translate [10, 0, 0], selection
```

### 이면각 설정

```python
# 이면각 직접 설정
set_dihedral \
  /prot/A/100/N, \
  /prot/A/100/CA, \
  /prot/A/100/C, \
  /prot/A/101/N, \
  180.0
```

### Undo/Redo

```python
undo
redo
```

---

## 무비 제작

### 프레임 설정

```python
# 프레임 수 설정
mset 1 x100              # 100 프레임

# 상태 기반
mset 1 -60               # 상태 1-60 사용
mset 1 x30 1 -60 60 x30  # 혼합
```

### 키프레임 애니메이션

```python
# 뷰 키프레임
mview store, 1           # 프레임 1에 현재 뷰 저장
# ... 뷰 변경 ...
mview store, 50          # 프레임 50에 뷰 저장
mview store, 100

# 보간
mview interpolate        # 중간 프레임 생성
```

### 씬 기반 애니메이션

```python
# Scene 전환 무비
scene S1, store
# ... 설정 변경 ...
scene S2, store

mset 1 x200
mview store, 1, scene=S1
mview store, 100, scene=S2
mview store, 200, scene=S1
mview interpolate
```

### 재생 제어

```python
mplay                    # 재생
mstop                    # 정지
mrewind                  # 처음으로
frame 50                 # 특정 프레임으로
forward                  # 앞으로
backward                 # 뒤로
```

### 내보내기

```python
# PNG 시퀀스
mpng frame_              # frame_0001.png, frame_0002.png, ...

# 레이트레이싱 적용
set ray_trace_frames, 1
mpng movie_frame_, ray=1
```

### FFmpeg로 비디오 변환

```bash
ffmpeg -framerate 30 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4
```

---

## 스크립팅

### 명령줄에서 스크립트 실행

```bash
pymol -c script.pml
pymol structure.pdb -d "run script.py"
```

### PML 스크립트

```python
# script.pml
load structure.pdb
hide everything
show cartoon
color spectrum
orient
ray 1920, 1080
png output.png
quit
```

### Python 스크립트

```python
# script.py
from pymol import cmd

cmd.load("structure.pdb")
cmd.hide("everything")
cmd.show("cartoon")
cmd.spectrum("count", "rainbow")
cmd.orient()
cmd.ray(1920, 1080)
cmd.png("output.png")
```

### Python API

```python
from pymol import cmd

# 명령 실행
cmd.do("load structure.pdb")

# 직접 호출
cmd.load("structure.pdb")
cmd.select("active", "resi 100-120")
cmd.show("sticks", "active")

# 정보 얻기
atoms = cmd.get_model("selection")
for atom in atoms.atom:
    print(atom.name, atom.resi, atom.coord)

# 좌표 얻기
coords = cmd.get_coords("selection")

# 거리 얻기
d = cmd.get_distance("/prot/A/10/CA", "/prot/A/20/CA")
```

### 사용자 정의 함수

```python
from pymol import cmd

def my_highlight(selection="all"):
    """잔기 하이라이트"""
    cmd.show("cartoon", "all")
    cmd.show("sticks", selection)
    cmd.color("yellow", selection)
    cmd.zoom(selection, 5)

# PyMOL에 등록
cmd.extend("highlight", my_highlight)

# 사용
# PyMOL> highlight resi 100-120
```

---

## 고급 기능

### 표면 계산

```python
# 분자 표면
show surface, protein

# 전기적 표면 (APBS 필요)
# 외부 도구로 .dx 파일 생성 후
load electrostatics.dx, emap
ramp_new eramp, emap, [-5, 0, 5]
set surface_color, eramp, protein
```

### 등고선/메쉬

```python
# 전자 밀도 맵
load map.ccp4, emap
isomesh mesh1, emap, 1.5          # 1.5σ 레벨
isomesh mesh2, emap, 3.0          # 3.0σ 레벨

# 스타일
color blue, mesh1
set mesh_width, 0.5
```

### 대칭 확장

```python
# 결정학적 대칭 mates
symexp prefix, object, selection, cutoff
symexp sym, protein, all, 5       # 5Å 이내 대칭 복사체
```

### 스테레오 뷰

```python
# 스테레오 모드
stereo on                 # 기본
stereo crosseye           # 교차시
stereo walleye           # 평행시
stereo quadbuffer        # 3D 안경 (하드웨어 필요)
stereo off
```

### 라벨링

```python
# 잔기 라벨
label n. CA, "%s-%s" % (resn, resi)

# 스타일
set label_color, black
set label_size, 14
set label_font_id, 7

# 위치 조정
set label_position, [0, 0, 2]
```

### 세션 저장

```python
# 전체 상태 저장
save session.pse

# 로드
load session.pse
```

---

## 명령어 레퍼런스

### 입출력

| 명령어 | 설명 |
|--------|------|
| `load` | 파일 로드 |
| `save` | 파일 저장 |
| `fetch` | PDB에서 다운로드 |
| `delete` | 객체 삭제 |
| `reinitialize` | 초기화 |

### 표현

| 명령어 | 설명 |
|--------|------|
| `show` | 표현 켜기 |
| `hide` | 표현 끄기 |
| `as` | 표현 교체 |
| `cartoon` | cartoon 스타일 |
| `label` | 라벨 표시 |

### 색상

| 명령어 | 설명 |
|--------|------|
| `color` | 색상 설정 |
| `spectrum` | 스펙트럼 색상 |
| `set_color` | 색상 정의 |
| `bg_color` | 배경색 |
| `recolor` | 다시 색칠 |

### 뷰

| 명령어 | 설명 |
|--------|------|
| `zoom` | 줌 |
| `center` | 중심 설정 |
| `orient` | 방향 정렬 |
| `origin` | 회전 중심 |
| `turn` | 회전 |
| `move` | 이동 |
| `clip` | 클리핑 |
| `view` | 뷰 저장/복원 |
| `scene` | 씬 관리 |

### 정렬

| 명령어 | 설명 |
|--------|------|
| `align` | 서열 기반 정렬 |
| `super` | 구조 기반 정렬 |
| `cealign` | CE 정렬 |
| `fit` | 좌표 피팅 |
| `pair_fit` | 원자쌍 피팅 |
| `rms` | RMSD 계산 |
| `rms_cur` | 현재 RMSD |

### 측정

| 명령어 | 설명 |
|--------|------|
| `distance` | 거리 측정 |
| `angle` | 각도 측정 |
| `dihedral` | 이면각 측정 |
| `get_area` | 표면적 |

### 편집

| 명령어 | 설명 |
|--------|------|
| `remove` | 원자 삭제 |
| `create` | 객체 생성 |
| `extract` | 추출 |
| `alter` | 속성 변경 |
| `h_add` | 수소 추가 |
| `rotate` | 회전 |
| `translate` | 이동 |
| `bond` / `unbond` | 결합 수정 |

### 렌더링

| 명령어 | 설명 |
|--------|------|
| `ray` | 레이트레이싱 |
| `draw` | 빠른 그리기 |
| `png` | PNG 저장 |
| `mpng` | 무비 PNG |

### 무비

| 명령어 | 설명 |
|--------|------|
| `mset` | 프레임 설정 |
| `mview` | 뷰 키프레임 |
| `mplay` / `mstop` | 재생/정지 |
| `frame` | 프레임 이동 |

### 설정

| 명령어 | 설명 |
|--------|------|
| `set` | 설정 변경 |
| `get` | 설정 조회 |
| `unset` | 설정 해제 |

---

## 유용한 팁

### 출판용 이미지

```python
# 고품질 설정
set ray_shadows, off
set antialias, 2
set ray_trace_mode, 1
bg_color white

# 렌더링
ray 4000, 3000
png figure.png, dpi=300
```

### 리간드-단백질 시각화

```python
# 기본 설정
hide everything
show cartoon, polymer
show sticks, organic
show sticks, byres organic around 4
color gray, polymer
util.cbay organic
distance hbonds, organic, byres organic around 4, mode=2
zoom organic, 8
```

### 여러 구조 비교

```python
# 로드
fetch 1abc 2def 3ghi

# 정렬
align 2def, 1abc
align 3ghi, 1abc

# 색상
color red, 1abc
color green, 2def
color blue, 3ghi

# 투명도
set cartoon_transparency, 0.5
```

---

## 참고 자료

- **공식 문서**: https://pymol.org/dokuwiki/
- **명령어 레퍼런스**: https://pymol.org/pymol-command-ref.html
- **PyMOL Wiki**: https://pymolwiki.org/
- **GitHub**: https://github.com/schrodinger/pymol-open-source
