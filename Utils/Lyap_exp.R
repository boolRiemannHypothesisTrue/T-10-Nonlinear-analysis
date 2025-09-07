# Установка и загрузка пакета (если не установлен)
# install.packages("tseriesChaos")
library(tseriesChaos)
setwd("C:/Users/MSI/Documents/GitHub/T-10-Nonlinear-analysis")
# --- Параметры ---
file_num <- 7
file_path <- paste0("T10_66131_Lpf", file_num, ".txt")

# Загрузка данных
time_series <- as.matrix(read.table(file_path))

time <- time_series[,1]
fp   <- time_series[,2]

# --- Analysis interval ---
time_min <- 400
time_max <- 900

mask <- (time >= time_min & time <= time_max)

time <- time[mask]
fp   <- fp[mask]

# --- Расчёт показателя Ляпунова ---
# Нужно выбрать embedding dimension (m) и задержку (d)
# Обычно подбирают через false.nearest и mutual
# Пусть fp = твой временной ряд
m   <- 5    # embedding dimension
d   <- 10   # задержка
s   <- 100  # число референсных точек
t   <- 20   # число шагов
ref <- 500  # длина временного окна для усреднения
k   <- 2    # количество ближайших соседей
eps <- 0.1  # радиус окрестности (подбирается!)

# вызов
lyap_result <- lyap_k(fp, m=3, d=2, s=200, t=40, ref=1700, k=2, eps=4)
# Посмотреть результат
print(lyap_result)
lyap(lyap_result, 5, 20)
recurr(fp, m=3, d=2, start.time=400, end.time=900)
stplot(fp, m=3, d=8, idt=1, mdt=250)

a <- d2(fp, m=4, d=2, t=4, eps.min=2)
