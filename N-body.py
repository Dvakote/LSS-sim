# -*- coding: utf-8 -*-
"""
Редактор Spyder

@author: Дмитрий Мелкозеров
"""

# v Подключаемые пакеты v
# ===========================================================================
import os
import math as m
import time
import numpy as np
import treecode.generate_system as generate
import treecode.create_tree_structure as cts
import treecode.tree_code as tc
import treecode.energy_and_momentum as EM
import treecode.make_graph as graph
import treecode.input_parametrs as parameters
# import threading
from joblib import Parallel, delayed
# ===========================================================================
# ^ Подключаемые пакеты ^
# v Используемые функции v
# ===========================================================================


def tree_root(Particles, Mass_center):
    # Функция, с которой начинается tree code
    if use_multiprocessing:
        # Весь объем системы разбивается на 8 кубов
        # каждый куб обсчитывается в своем процессе
        A0 = Parallel(n_jobs=workers, verbose=0)(
                delayed(tc.begin_tree)(Particles, Mass_center, i,
                                       n, eps_smooth)
                for i in range(1, 9))
        A = A0[0] + A0[1] + A0[2] + A0[3] + A0[4] + A0[5] + A0[6] + A0[7]
    else:
        A = np.zeros([np.size(Particles, 0), 4])
        # Весь объем системы разбивается на 8 кубов
        # каждый куб обсчитывается последовательно на одном ядре
        if not Mass_center[1, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 1, n, eps_smooth)
        if not Mass_center[2, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 2, n, eps_smooth)
        if not Mass_center[3, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 3, n, eps_smooth)
        if not Mass_center[4, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 4, n, eps_smooth)
        if not Mass_center[5, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 5, n, eps_smooth)
        if not Mass_center[6, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 6, n, eps_smooth)
        if not Mass_center[7, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 7, n, eps_smooth)
        if not Mass_center[8, 6] == 0:
            A += tc.begin_tree(Particles, Mass_center, 8, n, eps_smooth)
    return A


def tree_code_one_step(Y):
    # Метод leapfrog
    # рассчет изменения компонент скорости, где
    # Y[:, 3] - V_x, Y[:, 4] - V_y, Y[:, 5] - V_z,
    # Y[:, 7] - a_x, Y[:, 8] - a_y, Y[:, 9] - a_z;
    # по формуле V_{i + 1/2} = V_{i} + a_{i} * dt / 2
    Y[:, 3:6] += Y[:, 7:10] * time_step / 2
    # рассчет изменения координат, где
    # Y[:, 0] - x, Y[:, 1] - y, Y[:, 2] - z
    # по формуле X_{i + 1} = X_{i} + V_{i + 1/2} * dt
    Y[:, 0:3] += Y[:, 3:6] * time_step
    # Сортируем частицы в порядке возрастания номера ячейки
    Y = cts.distribution(Y, n, Distance)
    # Создаем древоподобную структуру "снизу вверх"
    Cells = cts.create_tree(Y, n, Distance)
    # Проходим древоподобную структуру "сверху вниз"
    # для получения ускорения и потенциала каждой частицы
    # Y[:, 10] - потенциал(phi) на i шаге
    Y[:, 7:11] = tree_root(Y, Cells)
    # Смотрим, есть ли материальные точки, вылетевшие
    # за рассматриваемый объем
    if Y[0, 11] < 0:
        # Считаем взаимодействие с улетевшими точками напрямую
        Y = tc.N_body_direct(Y, eps_smooth)
    Y[:, 7:11] *= G
    # рассчет изменения компонент скорости
    # по формуле V_{i + 1} = V_{i + 1/2} + a_{i + 1} * dt / 2
    Y[:, 3:6] += Y[:, 7:10] * time_step / 2
    return Y


def is_gravity_field_weak(Y):
    # Функция, выдающая ошибку, если гравитационное поле становится
    # слишком сильным для применения используемой модели
    global error
    global error_name
    # Считаем величину phi / c^2, где phi - гравитационный потенциал
    Array_phi = abs(Y[:, 10] / c_2)
    # если модуль phi / c^2 превышает определенное значение, то
    # гравитационное поле считаем сильным что выходит за границы
    # применимости используемой модели
    Array_phi = Array_phi >= 0.05
    if Array_phi.any():
        error = True
        error_name = 'Strong gravity field error'


def speed_limit(Y):
    # Функция, выдающая ошибку если скорость материальной
    # точки станет больше скорости света
    global error
    global error_name
    V = np.zeros([np.size(Y, 0), 3])
    V = np.multiply(Y[:, 3:6], Y[:, 3:6])
    V_2 = V.sum(axis=1) >= c_2
    if V_2.any():
        error = True
        error_name = 'FTL error'


if __name__ == "__main__":
    # v Константы v
    # =======================================================================
    # Гравитационная постоянная
    # G = 6.67408313 * m.pow(10, -11)  # м^3/(кг*с^2)
    G = 4.51811511 * m.pow(10, -15)  # кпк^3/(М_(Солнца)* (10^12 с)^2)
    # G = 4.51811511 * m.pow(10, -7)  # кпк^3/(М_(Млечного пути)* (10^15 с)^2)
    # Скорость света
    # c = 299792458 # м/с
    c = 9.7156188999  # кпк/(10^12 с)
# ===========================================================================
# ^ Константы ^
# v Параметры системы v
# ===========================================================================
# Прочие переменные (желательно не трогать)
    marker_size = 0.2  # 1
    c_2 = c * c
    error = False
    error_name = ''
    not_forbid_launch = True
    continue_input = True
    interrupted = False
    workers = os.cpu_count()

# Временной интервал
    # time_step = pow(10, 13)  # с
    time_step = 5.0  # 0.000025  # 10^12 с
    # time_step = 0.01  # 10^15 с

# Процентное распределение материи по типу
    d_e = 0.70  # Темная энергия
    d_m = 0.25  # Темная материя
    v_m = 0.05  # Видимая материя

# Параметр "сглаживания" гравитационного взаимодействия на близких дистанциях
    eps_smooth = 0.0  # кпк

# Параметры, которые нужны чаще всего (можно и нужно трогать)
# Количество ячеек по одной оси координат (для tree codes) в виде 2^(n)
    n = 4

# Минимальный размер ячейки по одной оси координат
    # Distance = 2 * 3.08567758 * pow(10, 22) # м
    Distance = 10 * m.pow(10, 3)  # кпк
    # Distance = 5 # Мпк

# Задаем первоначальный размер системы в единицах "Distance"
# для функции parameters_test
    i_test = 10
    j_test = 10
    k_test = 10
    indent_i = 0.0
    indent_j = 0.0
    indent_k = 0.0
    period = 0.0

# Параметры генерации эллипсоида в единицах (n * Distance / 2)
    a_inp = 1.0
    b_inp = 1.0
    c_inp = 1.0
# Начальные угловые скорости эллипсоида
    w_x = 0.0
    w_y = 0.0
    w_z = 0.0000005

# Средняя масса наблюдаемых объектов и их пекулярная скорость
    # m_avg = 1.98892 * pow(10, 41) # кг
    # v_avg = 0 #4 * pow(10, 5) / np.sqrt(3) # м/с
    m_avg = pow(10, 11)  # масс Солнц
    v_avg = 0.0  # 1.3 * pow(10, -2) / np.sqrt(3) # кпк/(10^12 c)
#     m_avg = 1 #масс Млечного пути
    # v_avg = 0 #1.3 * pow(10, -2) / np.sqrt(3) # Мпк/(10^15 c)

# Количество частиц
    N = 10000
# Число шагов
    Steps = 1
# Номера шагов, на которых требуется "сфотографировать положение всех
# материальных точек
    make_prelaunch_screenshot = False
    scr_step = []
# Тип сгенерированной системы (обязательно заполнить!)
    system_generation_type = 'ellipsoid'
# Использовать несколько процессов для вычислений
    use_multiprocessing = False
# Использовать данные, введенные вручную
    use_manual_input = False
# Использовать телеметрию
    use_telemetry = False
# Обратить время вспять
    inverse_time = False
# ===========================================================================
# ^ Параметры системы ^
# v Область с исполняемым кодом v
# ===========================================================================
    if use_manual_input:
        # Осуществляем ввод всех используемых параметров вручную
        inp_parameters = parameters.manual_input(c)
        print('Делать скриншоты системы?')
        print('y/n')
        input_variable = input()
        if (input_variable == 'y') or (input_variable == 'n'):
            if input_variable == 'y':
                make_prelaunch_screenshot = True
                enable_screenshots = True
                print('Наберите номера шагов, на которых нужно')
                print('сделать снимок системы')
                print('После того, как все нужные номера введены,')
                print('наберите "end" без кавычек, чтобы продолжить')
                while enable_screenshots:
                    input_var = input()
                    if input_var == 'end':
                        enable_screenshots = False
                    else:
                        try:
                            temp_var = int(input_var)
                            scr_step.append(temp_var)
                        except ValueError:
                            print('Номер шага может быть только целым числом')
            else:
                make_prelaunch_screenshot = False
                scr_step = []
        else:
            print('Введено недопустимое значение')
            print('Создание скриншотов отменено')
            make_prelaunch_screenshot = False
        print('Использовать телеметрию?')
        print('y/n')
        input_variable = input()
        if (input_variable == 'y') or (input_variable == 'n'):
            use_telemetry = input_variable == 'y'
        else:
            print('Введено недопустимое значение')
            print('Телеметрия не используется')
            use_telemetry = False
        print('Использовать многоядерность?')
        print('y/n')
        input_variable = input()
        if (input_variable == 'y') or (input_variable == 'n'):
            use_multiprocessing = input_variable == 'y'
        else:
            print('Введено недопустимое значение')
            print('Многоядерность не используется')
        print('Изменить знак у временного интервала?')
        print('y/n')
        input_variable = input()
        if (input_variable == 'y') or (input_variable == 'n'):
            inverse_time = input_variable == 'y'
        else:
            print('Введено недопустимое значение')
            inverse_time = False
        time_step = int(inp_parameters[18])
        Steps = int(inp_parameters[19])
        eps_smooth = int(inp_parameters[20])
        system_generation_type = str(inp_parameters[21])
    else:
        inp_parameters = [N, m_avg, v_avg,
                          n, Distance,
                          a_inp, b_inp, c_inp,
                          w_x, w_y, w_z,
                          i_test, j_test, k_test,
                          indent_i, indent_j, indent_k,
                          period]
    if (d_e >= 0) and (d_m >= 0) and (v_m > 0) \
            and (abs(1 - d_e - d_m - v_m) < 0.00000000001):
        m_avg = m_avg * (1 + (d_m / v_m))
    else:
        not_forbid_launch = False
        print('Недопустимое соотношение типов материи')
    if (time_step <= 0) or (Distance <= 0):
        not_forbid_launch = False
        print('Недопустимые параметры системы')
    if n > 0:
        n = int(m.pow(2, int(n)))
    else:
        not_forbid_launch = False
        print('Количество ячеек не может быть нулевым или отрицательным')
    if inverse_time:
        time_step *= -1
    try:
        try:
            try:
                # Генерируем или загружаем конфигурацию системы
                # в зваисимости от указанных/введенный параметров
                if (system_generation_type == 'cube')\
                    or (system_generation_type == 'random')\
                        or (system_generation_type == 'ellipsoid'):
                    X = generate.generate_system(system_generation_type,
                                                 inp_parameters)
                    if isinstance(X, str):
                        not_forbid_launch = False
                        print('Полуоси эллипсоида не могут быть нулевыми')
                    else:
                        np.savetxt('last config.txt', X)
                elif system_generation_type == 'last':
                    X = np.loadtxt('last config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif system_generation_type == 'debug':
                    X = np.loadtxt('error config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif system_generation_type == 'test':
                    X = np.loadtxt('test config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif system_generation_type == 'final':
                    X = np.loadtxt('final config.txt', dtype='float64')
                    N = np.size(X, 0)
                else:
                    not_forbid_launch = False
                    print('Выбранная конфигурация не может быть загружена')
            except IOError:
                not_forbid_launch = False
                print('Отсутствует необходимый файл конфигурации')
        except TypeError:
            not_forbid_launch = False
            print('Число материальных точек всегда должно быть целым')
    except ValueError:
        not_forbid_launch = False
        print('Неприемлимое число материальных точек')
    if not_forbid_launch:
        # Если загруженная конфигурация старая, то придаем ей
        # используемую в настоящий момент форму
        if np.size(X, 1) == 12:
            migration = np.zeros([np.size(X, 0), 2])
            X = np.hstack((X, migration))
            np.savetxt('last config.txt', X)
        # Определяем количество используемых процессов в зависимости
        # от используемого процессора
        if workers >= 8:
            workers = 8
        elif workers >= 4:
            workers = 4
        elif workers >= 2:
            workers = 2
        else:
            use_multiprocessing = False
        try:
            if make_prelaunch_screenshot:
                graph.screenshot(X, 'Шаг 0', marker_size)
            Energy = np.zeros([Steps, 6])
            start = time.time()
            for q in range(Steps):
                # Проверяем, не появились ли в системе сильные
                # гравитационные поля или скорости, равные или
                # превышающие световую
                speed_limit(X)
                is_gravity_field_weak(X)
                if error:
                    np.savetxt('error config.txt', X)
                    graph.screenshot(X, error_name, marker_size)
                    print(error_name + ' at step ' + str(q))
                    break
                # Выполняем один шаг по алгоритму tree code
                X = tree_code_one_step(X)
                # Записываем энергию системы в отдельный массив
                # 1) номер шага
                # 2) кинетическая энергия всей системы
                # 3) потенциальная энергия всей системы
                # 4) полная энергия всей системы
                # 5) максимальная разница в кинетической энергии
                # между исследуемым шагом и предыдущим
                # 6) максимальная разница в потенциальной энергии
                # между исследуемым шагом и предыдущим
                Energy[q] = [q,
                             EM.system_kinetic_energy(X),
                             EM.system_potential_energy(X),
                             EM.system_energy_Newton(X),
                             EM.max_dT(X),
                             EM.max_dU(X)]
                # записываем кинетическую и потенциальную энергию
                # каждой частицы для последующего использования
                X[:, 12] = EM.kinetic_energy_Newton(X)
                X[:, 13] = EM.potential_energy_Newton(X)
                if q in scr_step:
                    graph.screenshot(X, 'Шаг ' + str(q), marker_size)
            computing_time = time.time() - start
            print("Время выполнения", computing_time, "с")
            print(np.size(X, 0))
            if use_telemetry:
                EM.momentum_of_system(X)
                graph.full_telemetry(Energy, N)
        except KeyboardInterrupt:
            print('Работа программы прервана')
            EM.momentum_of_system(X)
            graph.full_telemetry(Energy, N)
        print('Сохранить финальную конфигурацию системы?')
        print('y?')
        input_variable = input()
        if input_variable == 'y':
            np.savetxt('final config.txt', X)
        else:
            print('Конфигурация не будет сохранена')
# ===========================================================================
# ^ Область с исполняемым кодом ^
