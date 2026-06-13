# 远程服务器与全球丰度分析入门指南

本文档面向第一次使用远程服务器做生信分析的使用者，目标是帮助你理解：

- 远程服务器、PyCharm、本地电脑之间是什么关系；
- 为什么要配置 conda 环境；
- 为什么全球 Candidatus Solibacter abundance 需要用 reads mapping 计算；
- 如何在服务器上跑通 3 个样本的 pilot mapping；
- 如何判断结果是否正常，以及常见错误怎么排查。

当前项目的远程服务器工作目录固定为：

```bash
/root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

如果后续服务器路径发生变化，需要把本文档中的路径替换成新的实际路径。

## 1. 远程服务器的基本运行逻辑

### 1.1 本地电脑、PyCharm 和远程服务器分别负责什么

可以把整个分析理解为三层：

```text
本地电脑
→ PyCharm
→ 远程服务器
```

本地电脑主要负责保存你的项目文件、脚本和文档。例如你本地的项目目录是：

```text
D:\pycharm\全球丰度分析1
```

PyCharm 是连接本地电脑和远程服务器的工具。你可以在 PyCharm 里查看文件、上传文件、打开远程 Terminal，并让远程服务器执行命令。

远程服务器是真正运行分析的地方。大型数据下载、MAG 解压、SRA reads 下载、CoverM mapping 都是在远程服务器上完成的，不是在本地电脑上完成的。

### 1.2 为什么要使用远程服务器

全球丰度分析涉及大量宏基因组数据，普通电脑通常会遇到几个问题：

- 磁盘不够；
- 下载速度慢；
- 内存不够；
- CPU 核心数少；
- 长时间运行容易中断。

远程服务器的优势是：

- 有更大的数据盘；
- 有更多 CPU 和内存；
- 可以持续运行较长时间；
- 适合下载和处理 SRA、MAG、宏基因组 reads 等大文件。

### 1.3 系统盘和数据盘的区别

服务器上通常有系统盘和数据盘。

系统盘用于安装系统和基础软件，不适合放大数据。你的服务器中数据盘是：

```bash
/root/autodl-tmp
```

因此本项目放在：

```bash
/root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

这样做的原因是避免把系统盘占满。

常用检查命令：

```bash
df -h
```

这个命令用于查看每个磁盘还剩多少空间。

## 2. PyCharm 中本地项目和远程项目的关系

你的本地 PyCharm 项目是：

```text
D:\pycharm\全球丰度分析1
```

你需要把本地的工作流文件夹上传到远程服务器：

```text
Result5_SMAG_global_abundance
```

上传后，远程服务器上的实际路径是：

```bash
/root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

在 PyCharm 里右键文件夹上传时，应选择：

```text
Deployment
→ Upload to ...
```

注意：

- 本地文件修改后，不会自动出现在服务器上，除非你上传或同步；
- 服务器上运行产生的大文件，也不会自动回到本地，除非你下载；
- 运行命令时，要确认 Terminal 连接的是远程服务器，而不是本地 Windows。

## 3. 项目目录结构说明

当前工作流目录中主要包含以下文件夹：

```text
Result5_SMAG_global_abundance
├── config
├── data
├── logs
├── results
├── scripts
├── README.md
├── READS_MAPPING_STEPS.md
├── REMOTE_PILOT_3_SAMPLES.md
└── REMOTE_SERVER_AND_GLOBAL_ABUNDANCE_GUIDE.md
```

各文件夹含义如下。

### 3.1 scripts

```text
scripts/
```

这里放自动化脚本。你后续主要运行这些脚本：

```bash
bash scripts/14_remote_pilot_preflight.sh
bash scripts/07_install_mapping_env.sh
bash scripts/08_download_smag_mag_archive.sh
bash scripts/09_extract_solibacter_mags_hpc.sh
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

### 3.2 config

```text
config/
```

这里放配置文件和样本列表。3 个 pilot 样本在：

```text
config/pilot_sra_accessions.txt
```

当前包含：

```text
SRR3984960
SRR3985396
SRR11833744
```

### 3.3 data/raw

```text
data/raw/
```

这里放下载的原始数据，例如：

- SMAG MAG 压缩包和解压文件；
- SRA reads；
- 原始 fastq 文件。

这是最占空间的目录。

### 3.4 data/processed

```text
data/processed/
```

这里放处理后的中间数据。例如从 SMAG MAG 中提取出来的 Candidatus Solibacter MAG reference set：

```text
data/processed/solibacter_mags
```

### 3.5 results

```text
results/
```

这里放结果文件。3 个 pilot 样本 mapping 后，最重要的结果是：

```text
results/sample_solibacter_abundance_pilot.tsv
```

这个表就是 3 个样本的 Candidatus Solibacter abundance 初版结果。

### 3.6 logs

```text
logs/
```

这里放运行日志。如果脚本报错，优先查看这里面的日志文件。

## 4. conda 环境是什么，为什么要激活

### 4.1 conda 环境的作用

生信分析需要很多命令行工具，例如：

- `coverm`
- `bwa`
- `samtools`
- `sra-tools`
- `seqkit`
- `pigz`

这些工具不是 Linux 系统自带的，需要安装。conda 环境的作用就是把这些工具集中安装到一个独立环境中，避免和系统环境冲突。

本项目使用的环境名是：

```bash
solibacter_mapping
```

### 4.2 激活 conda 环境

进入远程服务器后，先运行：

```bash
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
```

激活成功后，命令行前面通常会出现：

```text
(solibacter_mapping)
```

这表示你已经进入正确环境。

### 4.3 各工具的作用

`coverm`：计算 reads 回贴到 MAG 或基因后的 coverage 和 relative abundance。

`bwa`：把 reads 比对到参考基因组或 MAG。

`samtools`：处理比对后的 BAM 文件。

`sra-tools`：下载 NCBI SRA 数据，常用命令包括 `prefetch` 和 `fasterq-dump`。

`seqkit`：快速检查和处理 fasta/fastq 文件。

`pigz`：多线程压缩工具，比普通 `gzip` 更快。

## 5. 为什么要做 reads mapping

你想要得到的是：

```text
全球 Candidatus Solibacter abundance
```

这不能只靠一个 16S 序列检索来严格计算。原因是：

- 16S 检索只能说明某些样本中存在相似序列信号；
- MetaGraph 的 k-mer signal 更像是检出强度，不是真正的样本相对丰度；
- 要做正文级 abundance，最好回到 metagenomic reads，把 reads 比对到目标 MAG reference set。

因此更稳妥的路线是：

```text
SMAG MAG 数据库
→ 筛选 Candidatus Solibacter MAG
→ 构建 Solibacter reference set
→ 下载 SRA reads
→ reads mapping 到 Solibacter MAG
→ CoverM 计算 coverage / relative abundance
→ 汇总得到每个样本的 Candidatus Solibacter abundance
```

这个流程得到的结果比单纯 16S k-mer 检索更适合用于 Result 5。

## 6. 3 个样本 pilot mapping 的完整流程

### 6.1 进入项目目录

打开 PyCharm 的远程 Terminal 后，先运行：

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

确认当前位置：

```bash
pwd
```

应该显示：

```bash
/root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

### 6.2 激活 conda 环境

```bash
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
```

检查工具是否可用：

```bash
coverm --version
samtools --version
prefetch --version
fasterq-dump --version
seqkit version
```

如果这些命令能输出版本号，说明环境基本正常。

### 6.3 运行预检查

```bash
bash scripts/14_remote_pilot_preflight.sh
```

这个脚本会检查：

- 当前目录是否正确；
- 磁盘空间是否足够；
- pilot 样本列表是否存在；
- 必需的 metadata 和 manifest 是否存在；
- conda 和 mapping 工具是否安装；
- MAG 和 reads 是否已经下载。

### 6.4 下载 SMAG MAG 数据

```bash
bash scripts/08_download_smag_mag_archive.sh
```

这一步会下载 SMAG 的 MAG 压缩分卷并解压。它可能比较慢，因为 MAG 数据较大。

如果下载中断，可以重新运行同一个命令，通常会继续下载或跳过已完成文件：

```bash
bash scripts/08_download_smag_mag_archive.sh
```

检查下载情况：

```bash
ls -lh data/raw/smag_mag_parts
du -sh data/raw/smag_mag_parts
df -h /root/autodl-tmp
```

### 6.5 提取 Candidatus Solibacter MAG reference

MAG 下载完成后，运行：

```bash
bash scripts/09_extract_solibacter_mags_hpc.sh
```

检查是否提取成功：

```bash
ls data/processed/solibacter_mags | head
wc -l data/processed/solibacter_mag_ids.txt
```

如果能看到多个 MAG fasta 文件，说明 reference set 已经准备好。

### 6.6 下载 3 个 pilot 样本 reads

运行：

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
```

这一步会下载 3 个 SRA 样本：

```text
SRR3984960
SRR3985396
SRR11833744
```

检查 reads：

```bash
ls -lh data/raw/reads
du -sh data/raw/reads
```

### 6.7 运行 CoverM mapping

运行：

```bash
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

这一步会把 3 个样本的 reads 比对到 Candidatus Solibacter MAG reference set，并计算 abundance。

### 6.8 查看结果

运行：

```bash
cat results/sample_solibacter_abundance_pilot.tsv
```

如果这个文件存在，并且不是空文件，说明 3 个样本 pilot 跑通了。

## 7. 如何理解结果表

核心结果文件是：

```text
results/sample_solibacter_abundance_pilot.tsv
```

重点看以下几类列。

### 7.1 relative_abundance_sum

表示每个样本中 Candidatus Solibacter MAG 的相对丰度估计。

这是最接近文章中可视化 `Candidatus Solibacter abundance` 的指标。

### 7.2 trimmed_mean_sum

表示基于覆盖深度的稳健丰度指标。它会降低极端高覆盖区域对结果的影响。

### 7.3 covered_fraction_mean

表示目标 MAG 被 reads 覆盖的比例。

如果这个值非常低，说明虽然有少量 reads 命中，但证据可能不够强。

### 7.4 mapped_count_sum

表示比对到 Candidatus Solibacter MAG 的 reads 数量。

如果所有样本的 `mapped_count_sum` 都是 0，说明需要检查：

- Solibacter MAG 是否提取成功；
- reads 是否下载完整；
- mapping 参数是否过严；
- 这些样本中是否确实没有目标类群。

## 8. 常见错误和解决方法

### 8.1 No such file or directory

常见原因是当前不在正确目录。

先运行：

```bash
pwd
ls
```

如果没有看到 `scripts`、`config`、`results` 等目录，说明位置不对。

进入正确目录：

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

### 8.2 conda: command not found

说明 conda 没有被当前 shell 识别。

运行：

```bash
source /root/miniconda3/etc/profile.d/conda.sh
```

然后再激活环境：

```bash
conda activate solibacter_mapping
```

### 8.3 coverm: command not found

说明环境没有激活，或者工具没有安装成功。

先运行：

```bash
conda activate solibacter_mapping
```

再检查：

```bash
coverm --version
```

如果仍然找不到，重新安装环境：

```bash
bash scripts/07_install_mapping_env.sh
conda activate solibacter_mapping
```

### 8.4 下载没完成

MAG 或 reads 下载中断是正常情况，尤其是网络不稳定时。

重新运行原来的下载命令即可：

```bash
bash scripts/08_download_smag_mag_archive.sh
```

或：

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
```

### 8.5 磁盘不够

检查磁盘：

```bash
df -h
```

查看某个目录占用：

```bash
du -sh data/raw
du -sh data/raw/reads
du -sh data/raw/smag_mag_parts
```

如果空间不足，先不要跑 28 个样本或全量样本，只跑 3 个 pilot。

### 8.6 mapping 结果全是 0

如果结果表中 abundance 和 mapped count 全是 0，优先检查：

```bash
ls data/processed/solibacter_mags | head
ls -lh data/raw/reads
cat logs/extract_solibacter_mags.log
```

可能原因包括：

- Solibacter MAG reference 没有提取出来；
- reads 文件不完整；
- 样本本身不含目标类群；
- mapping identity 或 aligned percent 参数过严。

## 9. 每次上服务器后推荐的固定操作

每次重新打开服务器 Terminal 后，建议按这个顺序做：

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
pwd
df -h /root/autodl-tmp
bash scripts/14_remote_pilot_preflight.sh
```

这样可以快速确认：

- 当前目录正确；
- conda 环境正确；
- 磁盘空间够；
- 关键数据和工具状态正常。

## 10. MetaGraph signal 和 reads mapping abundance 的区别

前面已经用 MetaGraph 检索过 Ellin6076 16S sequence signal。

需要注意：

```text
MetaGraph 16S detection signal ≠ Candidatus Solibacter abundance
```

MetaGraph 的结果适合用于初步判断某个序列在全球样本中是否有检出信号，但它不是严格的丰度计算。

正文中如果要写：

```text
Candidatus Solibacter abundance
```

更推荐使用 reads mapping 结果，也就是本流程得到的：

```text
results/sample_solibacter_abundance_pilot.tsv
```

## 11. 当前阶段的目标

当前不要一开始就跑全量数据。推荐分三步：

```text
第一步：跑通 3 个 pilot 样本
第二步：扩展到 28 个候选 SRA 样本
第三步：再考虑全 SMAG 样本或更多全球样本
```

现在最重要的是第一步：

```bash
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

只要 3 个样本能成功生成：

```text
results/sample_solibacter_abundance_pilot.tsv
```

说明整个分析流程已经跑通。

## 12. 后续 pmtA abundance 怎么扩展

本文档重点讲 Candidatus Solibacter abundance。pmtA abundance 是后续扩展分析，需要另外构建 pmtA reference gene set。

后续 pmtA 分析大致流程是：

```text
Ellin6076 pmtA sequence
→ 搜索 SMAG MAG 中的 PmtA homologs
→ 构建 pmtA reference gene set
→ reads mapping 到 pmtA reference
→ 计算 pmtA abundance
→ 计算 pmtA abundance / Candidatus Solibacter abundance
```

在文章 Result 5 中，最理想的结果组合是：

- 全球 Candidatus Solibacter abundance 分布；
- 全球 pmtA abundance 分布；
- Candidatus Solibacter abundance 与 soil pH 的相关性；
- pmtA abundance 与 soil pH 的相关性；
- pmtA abundance / Candidatus Solibacter abundance 与 soil pH 的相关性。

但在当前阶段，先把 Candidatus Solibacter abundance 的 3 样本 pilot 跑通即可。
