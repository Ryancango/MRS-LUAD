# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/GOnewbee.R
# Role: GO enrichment used to annotate the candidate gene set and motivate the membrane-related signature terminology.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# 引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)GO
library(dplyr)
library("ggpubr")
library(ComplexHeatmap)

pvalueFilter=0.05        # p值过滤条件
p.adjustFilter=0.05      # 矫正后的p值过滤条件

# 定义图形的颜色
colorSel="p.adjust"
if(p.adjustFilter>0.05){ 
  colorSel="pvalue"
}

# 读取输入文件（修改为读取CSV，有表头）
rt=read.table("brown&down_genes.csv", header=TRUE, sep=",", check.names=F)     

# 提取交集基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        # 去除基因id为NA的基因
# gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

# GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<p.adjustFilter),]
# 输出显著富集的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

# 绘制 GO 圈图
ontology.col = c("#00CC33FF", "#FFC20AFF", "#CC33FFFF")  # 定义本体颜色
data = GO[order(GO$p.adjust), ]  # 按 p.adjust 排序
datasig = data[data$pvalue < 0.05, , drop = F]  # 筛选显著项
BP = datasig[datasig$ONTOLOGY == "BP", , drop = F]  # BP 类别
CC = datasig[datasig$ONTOLOGY == "CC", , drop = F]  # CC 类别
MF = datasig[datasig$ONTOLOGY == "MF", , drop = F]  # MF 类别
BP = head(BP, 6)  # 取前 6 个 BP
CC = head(CC, 6)  # 取前 6 个 CC
MF = head(MF, 6)  # 取前 6 个 MF
data = rbind(BP, CC, MF)  # 合并数据
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]  # 分配颜色

# 整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))  # 背景基因数
Gene = as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))  # 富集基因数
ratio = Gene / BgGene  # Rich Factor
logpvalue = -log(data$pvalue, 10)  # -log10(pvalue)
logpvalue.col = brewer.pal(n = 8, name = "Reds")  # 颜色梯度
f = colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = logpvalue.col)  # 颜色映射函数
BgGene.col = f(logpvalue)  # 应用颜色
df = data.frame(GO = data$ID, start = 1, end = max(BgGene))  # 扇区框架
bed2 = data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)  # pvalue 颜色条数据
bed3 = data.frame(GO = data$ID, start = 1, end = Gene, Gene = Gene)  # 基因数条数据（修正列名为 Gene 以清晰）
bed4 = data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)  # Rich Factor 数据
bed4$ratio = bed4$ratio / max(bed4$ratio) * 9.5  # 缩放 ratio

# 绘制圈图主体
pdf(file = "GO.circlize.brown.pdf", width = 10, height = 10)  # 输出 PDF
par(omi = c(0.1, 0.1, 0.1, 1.5))  # 设置边距
circos.par(track.margin = c(0.01, 0.01))  # circlize 参数
circos.genomicInitialize(df, plotType = "none")  # 初始化扇区
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {  # 绘制标签轨
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {  # 添加轴
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}

circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",  # pvalue 颜色条轨
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = value[, 2],
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], adj = 0, cex = 0.8, ...)
                    })

circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",  # 基因数条轨
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = '#BA55D3',
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], cex = 0.9, adj = 0, ...)
                    })

circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white", bg.col = "grey90",  # Rich Factor 轨
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for (j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ybottom = 0, ytop = value[, 1], col = value[, 2],
                                         border = NA, ...)
                    })

circos.clear()  # 清除 circlize

# 绘制中间图例
middle.legend = Legend(labels = c('Number of Genes', 'Number of Select', 'Rich Factor'),  # 修改为 Rich Factor（无0-1）
                       type = "points", pch = c(15, 15, 17), legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
                       title = "", nrow = 3, size = unit(3, "mm"))
circle_size = unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

# 绘制本体图例
main.legend = Legend(labels = c("Biological Process", "Cellular Component", "Molecular Function"),
                     type = "points", pch = 15, legend_gp = gpar(col = ontology.col), title_position = "topcenter",
                     title = "ONTOLOGY", nrow = 3, size = unit(3, "mm"), grid_height = unit(5, "mm"),
                     grid_width = unit(5, "mm"))

# 绘制 -log10(pvalue) 图例
logp.legend = Legend(labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
                     type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
                     title_position = "topcenter", grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
                     size = unit(3, "mm"))

lgd = packLegend(main.legend, logp.legend)  # 打包图例
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()  # 关闭设备

