# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/go气泡图三合一.R
# Role: Publication-style GO bubble plot for BP/CC/MF terms.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# ==============================================================================
# 0. 环境与数据准备
# ==============================================================================
# 假设 plot_data 已经在内存中 (包含 Description, Count, p.adjust, Category, GeneRatio_Num)
library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. 核心修复：重新排序 (从大到小)
# ==============================================================================
# 逻辑：先按 Category 分组，组内按 Count (数量) 升序排列
# 注意：ggplot 画图是从下往上画的，所以因子顺序要是 "小 -> 大"，画出来才是 "上大下小"
plot_data_sorted <- plot_data %>%
  group_by(Category) %>%
  arrange(Category, Count) %>% # 升序排列
  ungroup()

# 锁定因子顺序
plot_data_sorted$Description <- factor(plot_data_sorted$Description, 
                                       levels = unique(plot_data_sorted$Description))

# ==============================================================================
# 2. 绘制顶刊级气泡图 (修复图例 + 增强梯度)
# ==============================================================================
# 定义配色：使用 "Spectral" 色谱的变体，对比度极高
# 红色 = P值极小 (显著)，蓝色 = P值稍大
my_gradient <- c("#3288BD", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142")

p <- ggplot(plot_data_sorted, aes(x = GeneRatio_Num, y = Description)) +
  
  # --- 1. 画背景线 (Lollipop风格，增加引导性) ---
  geom_segment(aes(x = 0, xend = GeneRatio_Num, y = Description, yend = Description), 
               color = "grey90", size = 0.5) +
  
  # --- 2. 画气泡 ---
  # fill 映射 P值 (颜色)，size 映射 Count (大小)
  # shape=21 (实心圆带边框)，让气泡更有质感
  geom_point(aes(size = Count, fill = -log10(p.adjust)), 
             shape = 21,      # 实心圆
             color = "black", # 黑色边框
             stroke = 0.3,    # 边框细一点
             alpha = 0.9) +   # 不透明度
  
  # --- 3. 分面展示 ---
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  # --- 4. 颜色梯度 (增强版) ---
  scale_fill_gradientn(colours = my_gradient, 
                       name = "-log10(FDR)") + # 图例标题
  
  # --- 5. 大小范围 ---
  scale_size_continuous(range = c(3, 9), name = "Gene Count") +
  
  # --- 6. 主题美化 ---
  theme_bw() +
  theme(
    text = element_text(family = "sans", color = "black"),
    
    # 坐标轴
    axis.text.y = element_text(size = 11, color = "black"), # 通路名
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    
    # 分面标题 (放在右侧的灰色框)
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 10, face = "bold", color = "black", angle = 270),
    
    # 去掉杂乱网格，只留垂直辅助线
    panel.grid.major.y = element_line(color = "grey95", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 图例设置 (修复空白问题)
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  
  # --- 7. 修复图例显示 ---
  guides(
    # 颜色图例：做成连续的色条
    fill = guide_colorbar(order = 1, frame.colour = "black", ticks.colour = "black"),
    # 大小图例：强制显示灰色背景，解决“空白”问题
    size = guide_legend(order = 2, override.aes = list(fill = "grey70", color = "black"))
  ) +
  
  labs(x = "Gene Ratio", y = "")

# ==============================================================================
# 3. 保存
# ==============================================================================
# 自动计算高度
plot_height <- max(6, nrow(plot_data_sorted) * 0.35 + 2)

pdf("Fig1_GO_Bubble_Sorted_Gradient.pdf", width = 9.5, height = plot_height)
print(p)
dev.off()

print("✅ 修复完成！请查看 Fig1_GO_Bubble_Sorted_Gradient.pdf")