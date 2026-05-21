# Все библиотеки ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(corrplot)
library(car)
library(caret)
# Создание датасета в R -----
tbl <-"D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/Structured table tensie.csv"
tens <- read.csv2(tbl, header=TRUE, stringsAsFactors=FALSE, fileEncoding = "Windows-1250")
# Отделяем reference от Nb-B
ref_data <-tens[tens$Type == "Ref",]
nb_data <-tens[tens$Type == "Nb-B",]
# Тьюки------------------
stat <- tens %>%
  pivot_longer(cols = c(YS, Unknown, UTS, A),
               names_to = "Parameter",
               values_to = "Value") %>%
  group_by(Type, Size, Metrics, Parameter) %>%
  summarise(
    n = n(),
    Среднее = round(mean(Value), 3),
    Медиана = round(median(Value), 3),
    Стандартное_отклонение = round(sd(Value), 3),
    Дисперсия = round(var(Value), 3),
    Минимум = round(min(Value), 3),
    Максимум = round(max(Value), 3),
    Q1 = round(quantile(Value, 0.25), 3),
    Q3 = round(quantile(Value, 0.75), 3),
    IQR = round(IQR(Value), 3),
    Размах = round(diff(range(Value)), 3),
    Асимметрия = round(moments::skewness(Value), 3),
    Эксцесс = round(moments::kurtosis(Value), 3),
    .groups = "drop"
  )

print(stat)
# Гистограммы------
tens_long <- tens %>%
  select(Type, Metrics, Size, YS, Unknown, UTS, A) %>%
  pivot_longer(cols = c(YS, UTS, Unknown, A),
               names_to = "Parameter",
               values_to = "Value")

# Создаём комбинированную группу
tens_long$Group <- paste(tens_long$Type, tens_long$Metrics, sep = " - ")

calc_optimal_bins <- function(x, method = "sturges") {
  x_clean <- x[!is.na(x)]
  N <- length(x_clean)
  
  if (N < 2) return(5)  # минимальное значение
  
  bins <- switch(method,
                 sturges = round(log2(N) + 1),           # формула (14)
                 brooks = round(5 * log10(N)),            # формула (15)
                 sqrt = round(sqrt(N)),                   # формула (16)
                 cube_root = round(N^(1/3))               # формула (17)
  )
  
  # Ограничиваем разумными пределами
  bins <- max(3, min(bins, 30))
  return(bins)
}
# 4. Расчёт интервалов для каждой комбинации группы и параметра
tens_long <- tens_long %>%
  group_by(Group, Parameter) %>%
  mutate(
    N = n(),
    OptimalBins = calc_optimal_bins(Value, "sturges"),
    SturgesBins = calc_optimal_bins(Value, "sturges"),
    BrooksBins = calc_optimal_bins(Value, "brooks"),
    SqrtBins = calc_optimal_bins(Value, "sqrt"),
    CubeRootBins = calc_optimal_bins(Value, "cube_root")
  ) %>%
  ungroup()

# Просмотр результатов расчёта интервалов
cat("\n=== ОПТИМАЛЬНОЕ КОЛИЧЕСТВО ИНТЕРВАЛОВ ===\n")
bins_table <- tens_long %>%
  select(Group, Parameter, N, SturgesBins, BrooksBins, SqrtBins, CubeRootBins) %>%
  distinct()
print(bins_table)

G1 <- ggplot(tens_long, aes(x = Value, fill = Group)) +
  geom_histogram(bins = unique(tens_long$SturgesBins)[1], alpha = 0.7, position = "identity") +
  facet_grid(Group ~ Parameter, scales = "free") +
  scale_fill_manual(values = c("Ref - Flat" = "#3498db",
                               "Ref - Round" = "#2980b9",
                               "Nb-B - Flat" = "#e74c3c",
                               "Nb-B - Round" = "#c0392b")) +
  labs(title = "AM50 Magnesium Alloy: Mechanical Properties Distribution",
       subtitle = "Comparison: Reference vs Nb-B inoculation | Flat (5mm) vs Round (6mm)",
       x = "Value",
       y = "Frequency") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold"))
print(G1)
# Ящик с усами------  
tens$Group <- paste(tens$Type, tens$Size, tens$Metrics, sep = "_")
# Преобразуем факторы
tens$Type <- as.factor(tens$Type)
tens$Size <- as.factor(tens$Size)
tens$Metrics <- as.factor(tens$Metrics)
# Боксплот UTS по типу материала
BP1 <- ggplot(tens, aes(x = Type, y = UTS, fill = Size)) +
  geom_boxplot() +
  labs(x = "Тип материала", 
       y = "UTS (предел прочности, МПа)",
       title = "Сравнение предела прочности") +
  theme_minimal() +
  theme(legend.position = "none")
# Боксплот YS по типу материала
BP2 <- ggplot(tens, aes(x = Type, y = YS, fill = Size)) +
  geom_boxplot() +
  labs(x = "Тип материала", 
       y = "YS (предел текучести, МПа)",
       title = "Сравнение предела текучести") +
  theme_minimal()+
  theme(legend.position = "none")
# Боксплот Unknown по типу материала
BP3 <- ggplot(tens, aes(x = Type, y = Unknown, fill = Size)) +
  geom_boxplot() +
  labs(x = "Тип материала", 
       y = "Unknown",
       title = "Сравнение Unknown") +
  theme_minimal()+
  theme(legend.position = "none")
# Боксплот A по типу материала
BP4 <- ggplot(tens, aes(x = Type, y = A, fill = Size)) +
  geom_boxplot() +
  labs(x = "Тип материала", 
       y = "A (Относительное удлинение, 0,2%)",
       title = "Сравнение A")+
  theme_minimal() +
  theme(legend.position = "none")

(BP1|BP2)/(BP3|BP4)
# Проверка на нормальное распределение------
shap <- tens %>%
  group_by(Type, Size, Metrics) %>%
  summarise(across(
    c("YS", "Unknown", "UTS", "A"),
    ~ shapiro.test(.)$p.value))
print(shap)

# У всех параметров p-value> 0,05

# Диаграммы рассеяния------
# YS vs UTS
YU <- ggplot(tens, aes(x = YS, y = UTS, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS, МПа", 
       y = "UTS,МПа",
       title = "YS vs UTS") +
  theme_minimal()
# YS vs A
YA <- ggplot(tens, aes(x = YS, y = A, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS, МПа", 
       y = "A, 0,2%",
       title = "YS vs A") +
  theme_minimal()
# UTS vs  A
UA <- ggplot(tens, aes(x = UTS, y = A, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "UTS, МПа", 
       y = "A, 0,2%",
       title = "UTS vs A") +
  theme_minimal()
# YS vs Unknown
YUn <- ggplot(tens, aes(x = YS, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS, МПа", 
       y = "Unknown",
       title = "YS vs Unknown") +
  theme_minimal()
# UTS vs Unknown
UUn <- ggplot(tens, aes(x = UTS, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "UTS, МПа", 
       y = "Unknown",
       title = "UTS vs Unknown") +
  theme_minimal()

# A vs  Unknown
AUn <- ggplot(tens, aes(x = A, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "A, 0,2%", 
       y = "Unknown",
       title = "A vs Unknown") +
  theme_minimal()
YU/YA
UA/YUn
UUn/AUn
# Столбчатые диаграммы----
#UTS
stats_uts <- tens %>%
  group_by(Type, Size, Metrics) %>%
  summarise(
    Mean = mean(UTS),
    SD = sd(UTS),
    n = n(),
    SE = SD / sqrt(n),
    CI_lower = Mean - 2 * SE,
    CI_upper = Mean + 2 * SE,
    .groups = "drop"
  )
#UTS
MU <- ggplot(stats_uts, aes(x = interaction(Size, Metrics), y = Mean, fill = Type)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9),width = 0.3, color = "red") +
  labs(x = "Размер_Форма", 
       y = "Средний UTS (МПа)",
       title = "Средний предел прочности") +
  theme_minimal()
#YS
stats_uts <- tens %>%
  group_by(Type, Size, Metrics) %>%
  summarise(
    Mean = mean(YS),
    SD = sd(YS),
    n = n(),
    SE = SD / sqrt(n),
    CI_lower = Mean - 2 * SE,
    CI_upper = Mean + 2 * SE,
    .groups = "drop"
  )
MY <- ggplot(stats_uts, aes(x = interaction(Size, Metrics), y = Mean, fill = Type)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9),width = 0.3, color = "red") +
  labs(x = "Размер_Форма", 
       y = "Средний YS (МПа)",
       title = "Средний предел текучести") +
  theme_minimal()
#A
stats_uts <- tens %>%
  group_by(Type, Size, Metrics) %>%
  summarise(
    Mean = mean(A),
    SD = sd(A),
    n = n(),
    SE = SD / sqrt(n),
    CI_lower = Mean - 2 * SE,
    CI_upper = Mean + 2 * SE,
    .groups = "drop"
  )
MA <- ggplot(stats_uts, aes(x = interaction(Size, Metrics), y = Mean, fill = Type)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9),width = 0.3, color = "red") +
  labs(x = "Размер_Форма", 
       y = "Средний A (0,2% МПа)",
       title = "Среднее относительное удлинение") +
  theme_minimal()
#Unknown
stats_uts <- tens %>%
  group_by(Type, Size, Metrics) %>%
  summarise(
    Mean = mean(Unknown),
    SD = sd(Unknown),
    n = n(),
    SE = SD / sqrt(n),
    CI_lower = Mean - 2 * SE,
    CI_upper = Mean + 2 * SE,
    .groups = "drop"
  )
MUn <- ggplot(stats_uts, aes(x = interaction(Size, Metrics), y = Mean, fill = Type)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.9),width = 0.3, color = "red") +
  labs(x = "Размер_Форма", 
       y = "Средний Unknown",
       title = "Среднее Unknown") +
  theme_minimal()
(MU|MY)/(MA|MUn)
# Параметрические методы типа ANOVA ----
tens$Type <- as.factor(tens$Type)
tens$Size <- as.factor(tens$Size)
tens$Metrics <- as.factor(tens$Metrics)

params <- c("YS", "Unknown", "UTS", "A")
t_test_list <- list()

for (param in params) {
  t_test_list[[param]] <- t.test(as.formula(paste(param, "~ Type")), data = tens)
}

results_table <- data.frame(
  Parameter = params,
  Mean_NbB = round(sapply(t_test_list, function(x) x$estimate[1]), 2),
  Mean_Ref = round(sapply(t_test_list, function(x) x$estimate[2]), 2),
  Difference = round(sapply(t_test_list, function(x) x$estimate[1] - x$estimate[2]), 2),
  t_value = round(sapply(t_test_list, function(x) x$statistic), 3),
  p_value = sapply(t_test_list, function(x) format(x$p.value, scientific = TRUE, digits = 3)),
  Significant = ifelse(sapply(t_test_list, function(x) x$p.value < 0.001), "***",
                       ifelse(sapply(t_test_list, function(x) x$p.value < 0.01), "**",
                              ifelse(sapply(t_test_list, function(x) x$p.value < 0.05), "*", "ns")))
)

print(results_table)
# Корреляции (Пирсон) ----
cor_matrix <- cor(tens[, c("YS", "Unknown", "UTS", "A")])
corrplot(cor_matrix, method = "color", type = "full",
         diag = FALSE, addCoef.col = "white", number.cex = 1,
         title = "Корреляционная матрица")

# Линейная регрессия ----
#UTS
modelU <- lm(UTS ~ YS + A + Unknown, data = tens)
vif(modelU)
predicted <- predict(modelU)
actual <- tens$UTS
r2 <- R2(predicted, actual)
rmse <- RMSE(predicted, actual)
cat("R² =", round(r2, 4), "\n")
cat("RMSE =", round(rmse, 4), "\n")
#YS
modelY <- lm(YS ~ UTS+Unknown+A, data = tens)
vif(modelY)
predicted <- predict(modelY)
actual <- tens$YS
r2 <- R2(predicted, actual)
rmse <- RMSE(predicted, actual)
cat("R² =", round(r2, 4), "\n")
cat("RMSE =", round(rmse, 4), "\n")
#A
modelA <- lm(A ~ UTS+YS+Unknown, data = tens)
vif(modelA)
predicted <- predict(modelA)
actual <- tens$A
r2 <- R2(predicted, actual)
rmse <- RMSE(predicted, actual)
cat("R² =", round(r2, 4), "\n")
cat("RMSE =", round(rmse, 4), "\n")
#Unknown
modelUn <- lm(Unknown ~ UTS+A+YS, data = tens)
vif(modelUn)
predicted <- predict(modelUn)
actual <- tens$Unknown
r2 <- R2(predicted, actual)
rmse <- RMSE(predicted, actual)
cat("R² =", round(r2, 4), "\n")
cat("RMSE =", round(rmse, 4), "\n")
# Графика Моделей ----
#UTS
tens$predicted <- predict(modelU)
tens$residuals <- tens$UTS - tens$predicted
pred_U<-ggplot(tens, aes(x = predicted, y = UTS)) +
  geom_point(aes(color = abs(residuals)), size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Предсказанные значения UTS",
    y = "Реальные значения UTS",
    color = "|Ошибка|"
  ) +
  theme_minimal()
#A
tens$predicted <- predict(modelA)
tens$residuals <- tens$A - tens$predicted
pred_A<-ggplot(tens, aes(x = predicted, y = A)) +
  geom_point(aes(color = abs(residuals)), size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Предсказанные значения A",
    y = "Реальные значения A",
    color = "|Ошибка|"
  ) +
  theme_minimal()
#Un
tens$predicted <- predict(modelUn)
tens$residuals <- tens$Unknown - tens$predicted
pred_Un<-ggplot(tens, aes(x = predicted, y = Unknown)) +
  geom_point(aes(color = abs(residuals)), size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Предсказанные значения Unknown",
    y = "Реальные значения Unknown",
    color = "|Ошибка|"
  ) +
  theme_minimal()
#YS
tens$predicted <- predict(modelY)
tens$residuals <- tens$YS - tens$predicted
pred_Y<-ggplot(tens, aes(x = predicted, y = YS)) +
  geom_point(aes(color = abs(residuals)), size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Предсказанные значения YS",
    y = "Реальные значения YS",
    color = "|Ошибка|"
  ) +
  theme_minimal()

combo<-(pred_U|pred_Y)/(pred_A|pred_Un) +
  plot_annotation(
    title = "Визуализация RMSE"
  )
print(combo)

# Сохранение графиков ================================================

# 1. Гистограммы
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/histograms.png", 
       plot = G1, 
       width = 18, height = 12, 
       units = "cm", 
       dpi = 300)

# 2. Боксплоты
boxplots_combined <- (BP1 | BP2) / (BP3 | BP4)
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/boxplots.png", 
       plot = boxplots_combined, 
       width = 16, height = 14, 
       units = "cm", 
       dpi = 300)

# 3. Диаграммы рассеяния 
scatter_combined1 <- YU/YA
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/scatter_plots1.png", 
       plot = scatter_combined1, 
       width = 16, height = 12, 
       units = "cm", 
       dpi = 300)
scatter_combined2 <- UA/YUn
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/scatter_plots2.png", 
       plot = scatter_combined2, 
       width = 16, height = 12, 
       units = "cm", 
       dpi = 300)
scatter_combined3 <- UUn/AUn
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/scatter_plots3.png", 
       plot = scatter_combined3, 
       width = 16, height = 12, 
       units = "cm", 
       dpi = 300)

# 4. Корреляционная матрица
corrplot(cor_matrix, method = "color", type = "full",
         diag = FALSE, addCoef.col = "white", number.cex = 1,
         title = "Корреляционная матрица")
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/correlation_matrix.png", 
       plot=cormat,
       width = 10, height = 10, 
       units = "cm", 
       dpi = 300)

# 5. Визуализация RMSE (4 графика)
ggsave("D:/Dataset/Refinement of Mg alloys crystal structure_Brabant/rmse_plots.png", 
       plot = combo, 
       width = 18, height = 14, 
       units = "cm", 
       dpi = 300)
# Сохранение таблиц ----
write.csv2(stat, "D:/Dataset/Stats.csv", row.names = FALSE,  fileEncoding = "cp1251")
write.csv2(shap, "D:/Dataset/Shap.csv", row.names = FALSE,  fileEncoding = "cp1251")
write.csv2(results_table, "D:/Dataset/ttest.csv", row.names = FALSE,  fileEncoding = "cp1251")
