# check_gene · 单基因深度核查工具

对 rvtests-rvat 流程（`analysis/assoc_rvtest`）跑出的某个**基因**，把它的关联统计、变体明细、病例/对照携带情况、ToMMo 频率等汇总到一份报告里，方便对显著基因做人工核查。

- `check_gene.sh` —— 封装脚本：根据基因名 + impact 分层，自动定位流程输出文件并调用下面的 Python。
- `check_gene_detail.py` —— 核查主程序：解析 VCF / PLINK / pheno，计算每个变体、每个样本的细节。

---

## 用法

```bash
./check_gene.sh <Gene> [Impact] [Sample_Group]
```

| 参数 | 取值 | 默认 | 说明 |
|---|---|---|---|
| `Gene` | 基因符号，如 `STBD1` | （必填） | 与 rvtest `--geneFile`(refFlat) 里的基因名一致 |
| `Impact` | `high` / `moderate_high` / `low_moderate_high` | `moderate_high` | 变体影响分层（与流程的 `impactFilters` tag 对应；也兼容旧式 `impact_*` 输入） |
| `Sample_Group` | `both` / `case` / `control` | `both` | 输出哪一组的逐样本明细 |

### 示例

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_rvtest/check_gene

./check_gene.sh STBD1                       # 默认 moderate_high, both
./check_gene.sh STBD1 low_moderate_high     # 指定分层
./check_gene.sh STBD1 high case             # 只看 case 携带者
```

> 前提：对应的流程结果已经跑出来（`results/` 下有 `04.rvtest_run`、`05.post_process`、`03.info_filter` 等）。脚本会自动 `source activate cteph_geno_pro`。

---

## 它会读取哪些流程产物

脚本**按目录结构自动定位**（不写死长前缀，换 `inputPlinkPrefix` 也无需改脚本）：

| 内容 | 位置 |
|---|---|
| 原始 assoc（取 RANGE/NumVar） | `results/04.rvtest_run/<impact>/skato/*.SkatO.assoc` |
| CMC / SKAT-O / Zeggini 的 FDR 结果 | `results/05.post_process/<impact>/<method>/*.filtered.fdr.assoc` |
| 该分层的注释后 VCF | `results/03.info_filter/*.<impact>.vcf.gz` |
| 表型表 | `results/01.rvtest_prepare/*.pheno_rvt.tsv` |
| 基因模型 refFlat | `results/01.rvtest_prepare/refFlat.hg38.nochr.txt.gz` |
| 基因型矩阵（逐样本提取用） | `results/00.filter_genotype/*.gtqc.{bed,bim,fam}`（优先用 QC 过滤后的；若关闭了 `--filterGenotype`，回退到原始输入 PLINK） |
| ToMMo 60KJPN 等位基因频率 | `/LARGE0/.../ToMMo_60KJPN/tommo-60kjpn-...-autosome.norm.vcf.gz` |

> **重要**：基因型默认取 `00.filter_genotype` 的 **gtqc 过滤后** PLINK，所以报告里的样本数、携带计数、累积 MAC/MAF **都是去掉离群样本 + 去掉标记变体 + MAC≥2 之后的口径**，与正式分析一致。

---

## 输出

写到 `check_gene/output/<Gene>/` 下：

| 文件 | 内容 |
|---|---|
| `<Gene>.<impact>.summary.txt` | 基因级汇总：三个方法的 P/FDR、NumVar、累积 MAC/MAF（case vs ctrl）、**逐变体明细表**（impact/effect、case/ctrl 基因型计数、各自 AAF、ToMMo AAF、缺失率等） |
| `<Gene>.<impact>.case.sample_details.tsv` | 每个 **case** 样本携带了哪些变体、携带数、MAC |
| `<Gene>.<impact>.control.sample_details.tsv` | 每个 **control** 样本同上 |

终端还会彩色打印一份精简版（输入文件、基因统计框、变体表、逐样本携带者）。

### summary 关键字段

```
RVTest_NumVar / VCF_Hit_NumVar      rvtest 计入的变体数 / VCF 命中数
CMC_Pvalue / CMC_FDR                CMC burden 检验（rvtest 分数检验的 P/FDR）
SKAT-O_Pvalue / SKAT-O_FDR / Rho    SKAT-O 检验
Zeggini_Pvalue / Zeggini_FDR        Zeggini burden 检验（rvtest 分数检验的 P/FDR）
CMC_Effect_Beta/SE/OR/OR_95CI       CMC（携带指示，0/1）的效应量
Zeggini_Effect_Beta/SE/OR/OR_95CI   Zeggini（次要等位基因计数）的效应量
Cumulative_MAC_Gene (Case/Ctrl)     基因内累积次要等位基因计数
Cumulative_MAF_Case / _Ctrl         累积 MAF（= 累积 MAC / 2N），用于看病例富集
```

### 效应量（Beta / SE / OR / 95%CI）

rvtest 的 `cmc`/`zeggini` 是**分数检验，只给 P 值**。脚本额外用**逻辑回归**（协变量与流程一致：`sex + pc1_avg..pc10_avg`）算出两种折叠方式各自的效应量：

- **CMC**：burden = 是否携带任一罕见变异（指示变量 0/1）
- **Zeggini**：burden = 基因内累积次要等位基因数（计数）

模型 `logit P(case) = b0 + b·burden + 协变量`，`OR = exp(b)`，`95%CI = exp(b ± 1.96·SE)`。
不另报效应量的 p（用 rvtest 的 `CMC_Pvalue`/`Zeggini_Pvalue` 即可）。
> 已对拍：CMC 的 Beta/SE 与 rvtest 自带的 `cmcWald` 完全一致。两种 OR 在每个携带者只带 1 个等位基因时会相等，携带者带多个/纯合时才会分开。

---

## 依赖

- conda 环境 `cteph_geno_pro`（脚本会自动激活）
- `bcftools`、`plink2`（在该环境 PATH 中；plink2 路径写死在脚本里的 `PLINK2_PATH`）
- Python：`pandas`

---

## 常见调整

- **换基因型来源**：想对比"不过滤"的情形，把 `check_gene.sh` 里 `PLINK_PREFIX` 指到原始输入 PLINK，或先 `nextflow run rv_test.nf --filterGenotype false` 重跑。
- **换项目/路径**：改 `check_gene.sh` 顶部的 `ANALYSIS_DIR` / `TOMMO_VCF` / `INPUT_PLINK` 即可；文件名靠 glob 匹配，不依赖具体前缀。
- **新增检验方法**：流程里若再加方法（如 mb、fp），在 `check_gene_detail.py` 里仿照 `zeggini_stats` 加一段即可。
