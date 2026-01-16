# scRNA-seq Analysis Pipeline with Automatic PowerPoint Report Generation

## 개요

이 파이프라인은 scRNA-seq 데이터 분석을 자동화하고, 결과를 PowerPoint 프레젠테이션으로 자동 생성합니다.

### 주요 기능

✅ **자동화된 분석 파이프라인**
- Quality Control
- Normalization
- Batch Correction (Harmony/Scanorama)
- Clustering (Louvain/Leiden)
- Marker Gene Detection

✅ **자동 PowerPoint 리포트 생성**
- 모든 주요 figure 포함
- 통계 요약 테이블
- Marker gene 리스트
- 프레젠테이션 즉시 가능

✅ **완전한 결과 패키지**
- 처리된 AnnData 객체 (.h5ad)
- 고해상도 figure 파일 (PNG)
- Marker gene CSV 파일
- Metadata CSV 파일

---

## 설치

### 필수 패키지

```bash
pip install scanpy harmonypy python-igraph louvain python-pptx pandas numpy matplotlib seaborn
```

### 선택적 패키지

```bash
# Scanorama integration을 사용하려면
pip install scanorama

# Leiden clustering을 사용하려면
pip install leidenalg
```

---

## 사용법

### 1. 템플릿 사용 (새로운 분석)

템플릿 파일을 복사하고 설정을 수정합니다:

```bash
cp scrna_analysis_template.py my_analysis.py
```

`my_analysis.py` 파일을 열고 설정 섹션을 수정:

```python
# Analysis name (used for output files)
ANALYSIS_NAME = "My_Project_Name"

# Data paths
DATA_DIR = Path('/path/to/your/10x/data')
OUTPUT_DIR = Path('/path/to/output')

# Sample information
SAMPLES = {
    'Control_1': {'condition': 'Control', 'replicate': '1'},
    'Control_2': {'condition': 'Control', 'replicate': '2'},
    'Treatment_1': {'condition': 'Treatment', 'replicate': '1'},
    'Treatment_2': {'condition': 'Treatment', 'replicate': '2'},
}

# Analysis parameters (optional - defaults are reasonable)
PARAMS = {
    'min_genes': 200,
    'max_genes': 6000,
    'max_mt_percent': 20,
    'n_top_genes': 2000,
    'n_pcs': 50,
    'integration_method': 'harmony',  # 'harmony', 'scanorama', or None
    'clustering_method': 'louvain',   # 'louvain' or 'leiden'
    'resolutions': [0.5, 0.8, 1.0],
    'default_resolution': 0.8,
}
```

실행:

```bash
python my_analysis.py
```

### 2. 기존 결과에서 PowerPoint 생성

이미 분석이 완료된 데이터가 있는 경우:

```bash
python generate_ppt_report.py
```

또는 Python에서:

```python
from generate_ppt_report import create_scrna_report
from pathlib import Path

create_scrna_report(
    output_dir=Path('/path/to/output'),
    results_dir=Path('/path/to/results'),
    analysis_name="My Analysis"
)
```

---

## 출력 파일

분석 완료 후 다음 파일들이 생성됩니다:

```
output_directory/
├── Analysis_Name_Report.pptx          # 📊 PowerPoint 프레젠테이션
├── Analysis_Name_integrated.h5ad      # 💾 처리된 AnnData 객체
├── metadata.csv                        # 📋 Cell metadata
├── markers_cluster_0.csv               # 📊 Cluster 0 marker genes
├── markers_cluster_1.csv               # 📊 Cluster 1 marker genes
├── ...
└── figures/                            # 📁 모든 figure 파일
    ├── qc_metrics_before_filtering.png
    ├── highly_variable_genes.png
    ├── pca_variance_ratio.png
    ├── harmony_comparison.png          # (integration 사용 시)
    ├── umap_overview.png
    ├── cell_type_markers.png
    ├── marker_genes_per_cluster.png
    └── ...
```

---

## PowerPoint 리포트 내용

자동 생성되는 PowerPoint 프레젠테이션에는 다음이 포함됩니다:

### 슬라이드 구성

1. **Title Slide** - 분석 이름과 생성 날짜
2. **Analysis Overview** - 데이터셋 통계 요약
3. **Quality Control** - QC 메트릭 violin plots
4. **Highly Variable Genes** - 유전자 선택 결과
5. **Batch Correction** - Integration 전후 비교 (사용 시)
6. **PCA Variance** - 주성분 분산 비율
7. **UMAP Overview** - 전체 UMAP 시각화
8. **Condition Comparison** - 조건별 비교 분석
9. **Cell Distribution** - 조건별/클러스터별 세포 수 테이블
10. **Cluster Composition** - 클러스터 조성 분석
11. **Cell Type Markers** - 세포 타입 마커 발현
12. **Marker Genes** - 클러스터별 마커 유전자
13. **Top Markers per Cluster** - 주요 클러스터의 상위 마커 테이블
14. **Summary & Next Steps** - 요약 및 제안 사항

---

## 예제: JMV ALI 분석

### 실행 예제

```bash
# 전체 분석 실행 (PPT 자동 생성 포함)
python JMV_ALI_analysis.py

# 또는 기존 결과로 PPT만 생성
python generate_ppt_report.py
```

### 생성된 결과

```
JMV_ALI_results/
├── JMV_ALI_Organoid_scRNA-seq_Analysis_report.pptx  # 28 MB, 25+ slides
├── JMV_ALI_harmony_integrated.h5ad                   # 7.3 GB
├── metadata.csv                                       # 4.9 MB
├── markers_cluster_0.csv ~ markers_cluster_20.csv    # 21 files
└── figures/                                           # 15 PNG files
```

### 분석 결과

- **총 세포 수:** 35,245 cells (5개 조건)
- **클러스터:** 21 clusters (Louvain, resolution=0.8)
- **Integration:** Harmony (50 PCs)
- **실행 시간:** ~5-10분

---

## 커스터마이징

### Cell Type Markers 추가

조직 타입에 맞는 마커를 추가할 수 있습니다:

```python
CELL_TYPE_MARKERS = {
    # Airway epithelial cells
    'Basal cells': ['TP63', 'KRT5', 'KRT14'],
    'Ciliated cells': ['FOXJ1', 'PIFO', 'RSPH1'],
    'Club cells': ['SCGB1A1', 'SCGB3A2'],
    'Goblet cells': ['MUC5AC', 'MUC5B', 'TFF3'],

    # Immune cells
    'T cells': ['CD3D', 'CD3E', 'CD3G'],
    'B cells': ['CD79A', 'MS4A1'],
    'Macrophages': ['CD68', 'CD163', 'CSF1R'],
    'NK cells': ['NCAM1', 'KLRD1', 'NKG7'],
}
```

### PowerPoint 슬라이드 커스터마이징

`generate_ppt_report.py`의 `create_scrna_report` 함수를 수정하여:
- 슬라이드 추가/제거
- 레이아웃 변경
- 색상 스킴 수정
- 추가 통계 테이블 포함

---

## 통합 방법 비교

### Harmony (권장)
```python
PARAMS['integration_method'] = 'harmony'
```
- **장점:** 빠름, 효과적, 메모리 효율적
- **단점:** Python 환경에서 설치 필요
- **사용 시기:** 대부분의 경우 권장

### Scanorama
```python
PARAMS['integration_method'] = 'scanorama'
```
- **장점:** 매우 강력한 배치 효과 제거
- **단점:** 느림, 메모리 많이 사용
- **사용 시기:** 강한 배치 효과가 있을 때

### No Integration
```python
PARAMS['integration_method'] = None
```
- **사용 시기:** 단일 샘플 또는 배치 효과가 없을 때

---

## 클러스터링 방법 비교

### Louvain (기본값)
```python
PARAMS['clustering_method'] = 'louvain'
```
- **장점:** 빠름, 안정적
- **사용 시기:** 대부분의 경우

### Leiden
```python
PARAMS['clustering_method'] = 'leiden'
```
- **장점:** 더 정확한 커뮤니티 탐지
- **단점:** 약간 느림
- **사용 시기:** 세밀한 클러스터링이 필요할 때

---

## 문제 해결

### "ModuleNotFoundError: No module named 'scanpy'"

```bash
pip install scanpy
```

### "ModuleNotFoundError: No module named 'harmonypy'"

```bash
pip install harmonypy
```

### "ModuleNotFoundError: No module named 'igraph'"

```bash
pip install python-igraph louvain
```

### PowerPoint 생성 실패

```bash
pip install python-pptx
```

### 메모리 부족

- `n_pcs` 값을 줄입니다 (예: 30)
- `n_top_genes` 값을 줄입니다 (예: 1000)
- 샘플을 분할하여 처리합니다

---

## 성능 최적화

### 대용량 데이터셋 (>100K cells)

```python
PARAMS = {
    'n_top_genes': 2000,  # 유지
    'n_pcs': 30,          # 50 → 30으로 감소
    'n_neighbors': 10,    # 15 → 10으로 감소
}
```

### 빠른 테스트

```python
# 데이터 서브샘플링
adata = adata[adata.obs.sample(n=10000, random_state=42).index]
```

---

## 다음 단계

PowerPoint 리포트 생성 후:

1. **Cell Type Annotation**
   - 마커 유전자를 사용하여 클러스터에 세포 타입 할당
   - SingleR, CellTypist 등의 자동 annotation 도구 사용

2. **Differential Expression Analysis**
   - 조건 간 차등 발현 유전자 분석
   - Pseudobulk DESeq2 분석

3. **Pathway Enrichment**
   - GO, KEGG, Reactome 경로 분석
   - GSEA 분석

4. **Trajectory Analysis**
   - RNA velocity
   - Pseudotime analysis

5. **Integration with Other Data**
   - ATAC-seq, CITE-seq, Spatial transcriptomics

---

## 지원 및 문의

### 문서
- `JMV_ALI_ANALYSIS_SUMMARY.md` - 분석 결과 요약
- `README.md` - 결과 디렉토리 가이드

### 스크립트
- `scrna_analysis_template.py` - 범용 분석 템플릿
- `generate_ppt_report.py` - PPT 생성 스크립트
- `visualize_key_results.py` - 추가 시각화

---

**Last Updated:** 2026-01-15
**Version:** 1.0
