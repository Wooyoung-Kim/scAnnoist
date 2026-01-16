# Annotation Integration - 빠른 시작 (한국어)

## 개요

scRNA-agent는 이제 여러 annotation 방법을 통합하여 더 정확한 세포 타입 annotation을 생성할 수 있습니다.

**통합 가능한 방법:**
- 🔬 **CellTypist** - 딥러닝 기반, 빠르고 일관성 있음
- 📚 **Literature** - PubMed + CellMarker, 최신 마커 정보
- 🧬 **ScType** - 데이터베이스 기반, 학습 불필요

---

## 빠른 시작

### 1단계: CellTypist 실행

```python
from scrna_agent.tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)

# 데이터 로드
load_adata_for_annotation("data.h5ad")

# Tissue 정보 추가
add_tissue_metadata(tissue_map="tissue_map.csv", sample_col="sample")

# CellTypist annotation 실행
celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    output_dir="./outputs"
)
```

### 2단계: Literature-based Annotation 실행

```python
from scrna_agent.tools import annotate_with_literature_markers

# Literature 기반 annotation
annotate_with_literature_markers(
    tissue_type="Blood",
    cluster_key="leiden",
    output_col="literature_annotation"
)
```

### 3단계: 최적 Annotation 생성

```python
from scrna_agent.tools import create_optimal_annotation

# 두 방법을 통합하여 최적의 annotation 생성
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    confidence_threshold=0.5,
    output_col="optimal_annotation"
)

print(result)
```

✅ **완료!** `adata.obs['optimal_annotation']`에 최적의 annotation이 저장됩니다.

---

## 통합 전략

### 전략 1: Optimal Annotation (추천)

**특징:** Confidence와 일치도 기반 지능형 통합

**동작 방식:**
```
각 세포마다:
  CellTypist confidence >= threshold이면:
      → CellTypist 사용 (높은 신뢰도)
  그렇지 않으면:
      모든 방법의 예측 수집
      모든 방법이 일치하면:
          → Consensus 사용 (confidence = 0.9)
      그렇지 않으면:
          → 가중 투표 (weighted voting)
```

### 전략 2: Weighted Voting

**특징:** 각 방법에 가중치를 부여하여 투표

```python
from scrna_agent.tools import integrate_multiple_annotations

result = integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation", "sctype_annotation"],
    strategy="voting",
    weights={
        "celltypist_label": 1.0,
        "literature_annotation": 0.8,
        "sctype_annotation": 0.7
    },
    output_col="voting_result"
)
```

### 전략 3: Confidence-weighted

**특징:** Confidence 점수에 따라 자동으로 가중치 조정

```python
result = integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation"],
    strategy="confidence",
    confidence_cols={
        "celltypist_label": "celltypist_confidence"
    },
    output_col="confidence_integrated"
)
```

---

## 비교 및 분석

### 두 방법 비교

```python
from scrna_agent.tools import compare_celltypist_with_literature

comparison = compare_celltypist_with_literature(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation",
    cluster_key="leiden"
)

print(comparison)
```

**출력 예시:**
```
================================================================================
ANNOTATION COMPARISON REPORT
================================================================================

비교 방법: celltypist_label vs literature_annotation
전체 세포 수: 10,000

전체 일치도
-----------
일치율: 78.5%
불일치: 2,150 cells (21.5%)

클러스터별 일치도
----------------
Cluster 0: 92.3% 일치 (450/487 cells)
  CellTypist: B cells (95%)
  Literature: B cells (92%)

Cluster 1: 45.2% 일치 (123/272 cells)
  CellTypist: T cells (60%)
  Literature: NK cells (40%)
  ⚠ 높은 불일치 - 수동 검토 권장
...
```

---

## 완전한 워크플로우 예제

```python
from scrna_agent.tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    compare_celltypist_with_literature,
    create_optimal_annotation,
    plot_annotation_umap,
    save_annotations
)

# 1. 데이터 로드
print("[1/7] 데이터 로드...")
load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue_map="tissue_map.csv", sample_col="sample")

# 2. CellTypist annotation
print("[2/7] CellTypist annotation...")
celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    output_dir="./outputs"
)

# 3. Literature annotation
print("[3/7] Literature annotation...")
annotate_with_literature_markers(
    tissue_type="Blood",
    cluster_key="leiden",
    output_col="literature_annotation"
)

# 4. 비교
print("[4/7] 방법 비교...")
comparison = compare_celltypist_with_literature()
print(comparison)

# 5. 최적 annotation 생성
print("[5/7] 최적 annotation 생성...")
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    output_col="optimal_annotation"
)
print(result)

# 6. 시각화
print("[6/7] UMAP 시각화...")
plot_annotation_umap(
    annotation_col="optimal_annotation",
    output_file="./outputs/umap_optimal.png"
)

# 7. 저장
print("[7/7] 결과 저장...")
save_annotations(
    annotation_cols=[
        "celltypist_label",
        "literature_annotation",
        "optimal_annotation"
    ],
    output_file="./outputs/all_annotations.csv"
)

print("\n✓ 완료!")
print("  결과: ./outputs/all_annotations.csv")
print("  UMAP: ./outputs/umap_optimal.png")
```

---

## Python API

### 주요 함수

#### `create_optimal_annotation()`
여러 방법을 지능적으로 통합하여 최적의 annotation 생성

**파라미터:**
- `celltypist_col`: CellTypist annotation 컬럼 (필수)
- `celltypist_conf_col`: CellTypist confidence 컬럼
- `literature_col`: Literature annotation 컬럼 (선택)
- `sctype_col`: ScType annotation 컬럼 (선택)
- `confidence_threshold`: CellTypist를 직접 사용할 임계값 (기본: 0.5)
- `output_col`: 출력 컬럼 이름

**생성되는 컬럼:**
- `{output_col}`: 최적 annotation
- `{output_col}_confidence`: Confidence 점수
- `{output_col}_source`: 사용된 방법

#### `integrate_multiple_annotations()`
여러 방법을 투표 또는 confidence로 통합

**파라미터:**
- `methods`: 컬럼 이름 리스트
- `strategy`: "voting" 또는 "confidence"
- `weights`: 가중치 딕셔너리 (선택)
- `confidence_cols`: Confidence 컬럼 매핑 (선택)
- `output_col`: 출력 컬럼 이름

#### `compare_celltypist_with_literature()`
CellTypist와 Literature annotation 비교

**반환값:** 상세 비교 리포트 (문자열)

---

## 예제 파일

자세한 예제는 다음 파일 참조:
- **Python 예제**: `examples/annotation_integration_example.py`
- **전체 가이드**: `docs/ANNOTATION_INTEGRATION_GUIDE.md`

---

## 사용 시나리오

### 시나리오 1: CellTypist가 불확실한 경우

```python
# CellTypist로 먼저 annotation
celltypist_annotate_tissue_aware(output_dir="./outputs")

# Confidence가 낮은 세포 확인
import scanpy as sc
adata = sc.read_h5ad("annotated.h5ad")
low_conf_cells = adata.obs[adata.obs['celltypist_confidence'] < 0.5]
print(f"낮은 confidence: {len(low_conf_cells)} cells")

# Literature로 보완
annotate_with_literature_markers(tissue_type="Blood")

# 최적 annotation (낮은 confidence는 literature 활용)
create_optimal_annotation(
    confidence_threshold=0.5,
    output_col="optimal_annotation"
)
```

### 시나리오 2: 희귀 세포 타입

```python
# 세 가지 방법 모두 실행
celltypist_annotate_tissue_aware(output_dir="./outputs")
annotate_with_literature_markers(tissue_type="Blood", top_n_genes=100)
sctype_annotate(tissue="Blood")

# 세 방법 통합 (희귀 세포 타입 포착)
integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation", "sctype_annotation"],
    strategy="voting",
    output_col="comprehensive_annotation"
)
```

### 시나리오 3: 조직별 특수 세포

```python
# Tissue-specific literature search
annotate_with_literature_markers(
    tissue_type="Liver",  # 간 특이적
    top_n_genes=50,
    output_col="liver_literature"
)

# CellTypist와 통합
create_optimal_annotation(
    celltypist_col="celltypist_label",
    literature_col="liver_literature",
    confidence_threshold=0.3,  # 더 많이 literature 활용
    output_col="liver_optimal"
)
```

---

## 문제 해결

### 문제: 방법 간 일치도가 낮음 (<60%)

**원인:**
- 다른 조직/샘플 타입 사용
- 세포 타입 세분화 수준 차이
- 잘못된 tissue model

**해결:**
```python
# Tissue model 확인
from scrna_agent.tools.model_management import get_models_for_adata
models = get_models_for_adata(adata)
print(models)

# 더 관대한 통합
create_optimal_annotation(
    confidence_threshold=0.3,  # 낮은 임계값
    output_col="permissive_annotation"
)
```

### 문제: "Unknown" 세포가 많음 (>20%)

**해결:**
```python
# 더 많은 마커 사용
annotate_with_literature_markers(
    top_n_genes=100,  # 기본값 50에서 증가
    min_confidence=0.2  # 기본값 0.3에서 감소
)

# ScType 추가
sctype_annotate(tissue="Blood")

# 세 방법 통합
create_optimal_annotation(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation",
    sctype_col="sctype_annotation"
)
```

---

## 요약

**설치:** scRNA-agent에 이미 포함됨

**기본 사용법:**
```python
from scrna_agent.tools import (
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    create_optimal_annotation
)

# 두 방법 실행
celltypist_annotate_tissue_aware(output_dir="./outputs")
annotate_with_literature_markers(tissue_type="Blood")

# 통합
create_optimal_annotation(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation"
)
```

**주요 장점:**
- ✅ Consensus를 통한 높은 정확도
- ✅ Edge case 처리 개선
- ✅ Confidence 점수 제공
- ✅ 충돌 감지 및 해결

---

**이제 최적의 annotation을 만들 수 있습니다!** 🎯

## 더 알아보기

- **전체 가이드 (영문)**: `docs/ANNOTATION_INTEGRATION_GUIDE.md`
- **Python 예제**: `examples/annotation_integration_example.py`
- **Tissue-aware CellTypist**: `docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md`
- **PPT 생성**: `docs/PPT_QUICK_START.md`
