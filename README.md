# scAnnoist

> **sc** (single-cell) + **Anno** (annotation) + **ist** (specialist)

AI-powered single-cell RNA-seq analysis pipeline with automated cell type annotation and reporting.

## 🎯 주요 기능

scAnnoist는 scRNA-seq 데이터 분석을 완전 자동화하고, PowerPoint 리포트를 자동 생성하는 multi-agent 시스템입니다.

### ✨ 핵심 특징

1. **완전 자동화된 분석**
   - Quality Control & Filtering
   - Normalization & Feature Selection
   - Batch Correction (Harmony/Scanorama)
   - Clustering (Louvain/Leiden)
   - Marker Gene Detection

2. **종합적인 Cell Type Annotation** 🆕
   - **CellTypist** 자동 annotation
   - **Marker 기반** 수동 annotation
   - **MCP-based 문헌 검증** (PubMed 통합) ✨ NEW
   - 두 방법의 **비교 분석**
   - **최종 결정 근거** 상세 문서화

3. **자동 PowerPoint 리포트 생성**
   - 모든 주요 figure 포함
   - Annotation 비교 및 근거
   - 통계 요약 테이블
   - Marker gene 리스트
   - 즉시 프레젠테이션 가능

4. **완전한 결과 패키지**
   - Annotated AnnData (.h5ad)
   - 고해상도 figures (PNG)
   - Annotation report (CSV, MD)
   - Marker genes (CSV)
   - PowerPoint presentation

---

## 📦 설치

```bash
pip install scanpy harmonypy python-igraph louvain python-pptx celltypist pandas numpy matplotlib seaborn
```

---

## 🚀 빠른 시작

### 1단계: 기본 분석
```bash
python JMV_ALI_analysis.py
```

### 2단계: Cell Type Annotation
```bash
python comprehensive_cell_annotation.py
```

### 3단계: PowerPoint 생성
```bash
python generate_ppt_report.py
```

---

## 📊 출력 결과

- **PowerPoint 리포트** (25+ slides, ~28 MB)
- **Annotated Data** (.h5ad format)
- **Annotation 문서** (ANNOTATION_REPORT.md)
- **모든 Figures** (고해상도 PNG)
- **통계 테이블** (CSV)

---

## 🔬 MCP 통합 (NEW)

### PubMed 문헌 기반 Annotation 검증

이 파이프라인은 MCP (Model Context Protocol)를 통해 PubMed 문헌 데이터베이스에 접근하여 cell type annotation을 검증합니다.

**주요 기능:**
- 📚 자동 문헌 검색으로 marker 유전자 검증
- ✅ 논문 근거 기반 annotation 신뢰도 평가
- 🔍 관련 마커 자동 발견
- 📄 전문(full-text) 논문 분석

**사용 예시:**
```python
from scrna_agent.tools import mcp_validate_cell_type_annotation

# Annotation 검증
validation = mcp_validate_cell_type_annotation(
    cell_type="CD8 T cells",
    observed_markers=["CD8A", "CD8B", "GZMB"],
    species="human",
    confidence_threshold=0.7
)
```

**자세한 사용법:** `docs/MCP_INTEGRATION_GUIDE.md` 참조

---

자세한 내용은 `ANALYSIS_GUIDE.md`를 참조하세요.
