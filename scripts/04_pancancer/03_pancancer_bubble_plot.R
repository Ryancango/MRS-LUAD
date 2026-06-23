# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/xena/泛癌分析气泡图.R
# Role: Generate the pan-cancer survival-association bubble plot.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

library(ggplot2)
library(dplyr)
library(scales)

# ==========================================================
# 1. 数据准备与标签逻辑构建
# ==========================================================

# A. 基础清洗
final_data <- plot_df %>%
  mutate(
    logP = -log10(P),
    logP_capped = ifelse(logP > 8, 8, logP),
    HR_capped = ifelse(HR > 3, 3, ifelse(HR < 1, 1, HR)),
    stars = case_when(
      P < 0.001 ~ "***",
      P < 0.01 ~ "**",
      P < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Type = factor(Type, levels = c("OS", "DSS", "DFI", "PFI"))
  )

# B. 构建 "癌症 -> 标签" 的字典
# 确保每个癌症对应的 N 是唯一的
cancer_labels_df <- final_data %>%
  select(Cancer, N) %>%
  distinct() %>%
  mutate(
    # 核心步骤：生成 "LGG (N=533)" 格式的字符串
    Axis_Label = paste0(Cancer, " (N=", N, ")")
  )

# C. 智能排序 (按 HR 从高到低)
risk_order_code <- final_data %>%
  group_by(Cancer) %>%
  summarise(mean_HR = mean(HR, na.rm = TRUE)) %>%
  arrange(desc(mean_HR)) %>%
  pull(Cancer)

# D. 将排序逻辑应用到新标签上
# 这一步很关键：找到 risk_order_code 对应的 Axis_Label
ordered_labels <- cancer_labels_df$Axis_Label[match(risk_order_code, cancer_labels_df$Cancer)]

# E. 将数据中的标签设为有序因子
final_data <- final_data %>%
  left_join(cancer_labels_df, by = c("Cancer", "N")) %>%
  mutate(Axis_Label = factor(Axis_Label, levels = ordered_labels))

# ==========================================================
# 2. 绘图代码 (X轴改为 Axis_Label)
# ==========================================================

p_final_N <- ggplot(final_data, aes(x = Axis_Label, y = Method)) +
  
  # --- 气泡层 ---
  geom_point(aes(fill = HR_capped, size = logP_capped), 
             shape = 21, color = "black", stroke = 0.25) +
  
  # --- 星号层 ---
  geom_text(aes(label = stars), vjust = 0.76, size = 2.5, fontface = "bold") +
  
  # --- 配色 (白-金-红-深红) ---
  scale_fill_gradientn(
    colors = c("#F7F7F7", "#FFEDA0", "#FEB24C", "#F03B20", "#BD0026"),
    values = rescale(c(1, 1.2, 1.5, 2, 3)), 
    name = "Hazard Ratio",
    limits = c(1, 3), 
    oob = scales::squish
  ) +
  
  # --- 大小图例 ---
  scale_size(
    range = c(2, 8), 
    name = "Significance\n(-log10 P)",
    breaks = c(1.3, 3, 5),
    labels = c("p<0.05", "p<0.001", "p<1e-5")
  ) +
  
  # --- 分面 ---
  facet_grid(Type ~ ., switch = "y") +
  
  # --- 主题精修 ---
  theme_minimal(base_family = "sans") +
  theme(
    # X轴文字：旋转90度，字体稍微调小一点以防拥挤
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                               color = "black", size = 9, face = "plain"), # 既然有数字，用plain可能比bold更清晰
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.title = element_blank(),
    
    # 网格与边框
    panel.grid.major = element_line(color = "grey92", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    strip.background = element_rect(fill = "grey95", color = "black", size = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    
    # 图例与布局
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  labs(
    title = "Landscape of Prognostic Risk (Sorted by HR Magnitude)",
    subtitle = "Sample sizes (N) are shown in parentheses",
    caption = "Source: TCGA Pan-Cancer Data"
  )

# ==========================================================
# 3. 输出与保存
# ==========================================================

print(p_final_N)

# 建议宽度设大一点 (width=16)，因为现在的标签变长了
ggsave("PanCancer_Bubble_With_N_Sorted.pdf", p_final_N, width = 16, height = 7.5)
cat("✅ 最终图表已保存: PanCancer_Bubble_With_N_Sorted.pdf\n")