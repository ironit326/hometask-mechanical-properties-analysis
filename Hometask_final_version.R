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
# Отделяем по типу образцов
rr <- ref_data[ref_data$Metrics == "Round", ]
rf <- ref_data[ref_data$Metrics == "Flat", ]
nr <- nb_data[nb_data$Metrics == "Round", ]
nf <- nb_data[nb_data$Metrics == "Flat", ]

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
       title = "Сравнение Unknown по типу материала") +
  theme_minimal()+
  theme(legend.position = "none")
# Боксплот A по типу материала
BP4 <- ggplot(tens, aes(x = Type, y = A, fill = Size)) +
  geom_boxplot() +
  labs(x = "Тип материала", 
       y = "A (Относительное удлинение, 0,2%)",
       title = "Сравнение A по типу материала")+
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
YU <- ggplot(df, aes(x = YS, y = UTS, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (предел текучести, МПа)", 
       y = "UTS (предел прочности, МПа)",
       title = "YS vs UTS") +
  theme_minimal()
theme(legend.position = "none")

# YS vs A
YA <- ggplot(df, aes(x = YS, y = A, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (предел текучести, МПа)", 
       y = "A (Относительное удлинение, 0,2%)",
       title = "YS vs A") +
  theme_minimal()
theme(legend.position = "none")
# UTS vs  A
UA <- ggplot(df, aes(x = UTS, y = A, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (UTS (предел прочности, МПа)", 
       y = "A (Относительное удлинение, 0,2%)",
       title = "UTS vs A") +
  theme_minimal()
theme(legend.position = "none")
# YS vs Unknown
YUn <- ggplot(df, aes(x = YS, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (предел текучести, МПа)", 
       y = "Unknown",
       title = "YS vs Unknown") +
  theme_minimal()
theme(legend.position = "none")
# UTS vs Unknown
UUn <- ggplot(df, aes(x = UTS, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (UTS (предел прочности, МПа)", 
       y = "A (Относительное удлинение, 0,2%",
       title = "UTS vs Unknown") +
  theme_minimal()
theme(legend.position = "none")

# A vs  Unknown
AUn <- ggplot(df, aes(x = A, y = Unknown, shape = Size, color = Type )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
  labs(x = "YS (UTS (предел прочности, МПа)", 
       y = "Unknown",
       title = "A vs Unknown") +
  theme_minimal()
theme(legend.position = "none")
(YU|YA)/UA
(YUn|UUn)/AUn
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

anova_results <- tens %>%
  group_by(Size, Metrics) %>%
  summarise(
    # ANOVA для YS
    YS_anova = list(aov(YS ~ Type, data = pick(everything()))),
    # ANOVA для UTS
    UTS_anova = list(aov(UTS ~ Type, data = pick(everything()))),
    # ANOVA для A
    A_anova = list(aov(A ~ Type, data = pick(everything()))),
    # ANOVA для Unknown (убрали запятую в конце!)
    Unknown_anova = list(aov(Unknown ~ Type, data = pick(everything()))),
    .groups = "drop"
  )

anova_summary <- anova_results %>%
  mutate(
    YS_pvalue = sapply(YS_anova, function(x) summary(x)[[1]][["Pr(>F)"]][1]),
    UTS_pvalue = sapply(UTS_anova, function(x) summary(x)[[1]][["Pr(>F)"]][1]),
    A_pvalue = sapply(A_anova, function(x) summary(x)[[1]][["Pr(>F)"]][1]),
    Unknown_pvalue = sapply(Unknown_anova, function(x) summary(x)[[1]][["Pr(>F)"]][1])
  ) %>%
  select(Size, Metrics, YS_pvalue, UTS_pvalue, A_pvalue, Unknown_pvalue)
print(anova_summary)
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
modelA <- lm(A ~ UTS+Unknown+YS, data = tens)
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
