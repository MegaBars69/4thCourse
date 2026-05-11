#pragma once

// Запуск таблицы времени стабилизации (пункт 1 отчёта)
void task3_steady_time_table();

// Генерация данных для визуализации (пункт 2)
void task3_visualization(double w, double rho_in, int ramp_steps = 0);

// Таблица зависимости времени стабилизации от параметров (пункт 3)
void task3_convergence_table();