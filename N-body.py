# -*- coding: utf-8 -*-
"""
Редактор Spyder

@author: Дмитрий Мелкозеров
"""

import os
import math as m
import time
import numpy as np
from joblib import Parallel, delayed
import treecode.generate_system as generate
import treecode.create_tree_structure as cts
import treecode.tree_code as tc
import treecode.energy_and_momentum as EM
import treecode.make_graph as graph
import treecode.input_parametrs as parameters
import treecode.telemetry as tlm
import treecode.fractal as fractal
# import threading


def tree_root(particles, mass_center):
    # Функция, с которой начинается tree code
    if USE_MULTIPROCESSING:
        # Весь объем системы разбивается на 8 кубов
        # каждый куб обсчитывается в своем процессе
        a_parallel = Parallel(n_jobs=WORKERS, verbose=0)(
            delayed(tc.begin_tree)(particles, mass_center, i+1,
                                   SYSTEM_PARAMETERS)
            for i in range(0, 8))
        accel = a_parallel[0] + a_parallel[1] + a_parallel[2] + a_parallel[3] \
            + a_parallel[4] + a_parallel[5] + a_parallel[6] + a_parallel[7]
    else:
        accel = np.zeros([np.size(particles, 0), 4])
        # Весь объем системы разбивается на 8 кубов
        # каждый куб обсчитывается последовательно на одном ядре
        if not mass_center[1, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 1,
                                   SYSTEM_PARAMETERS)
        if not mass_center[2, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 2,
                                   SYSTEM_PARAMETERS)
        if not mass_center[3, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 3,
                                   SYSTEM_PARAMETERS)
        if not mass_center[4, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 4,
                                   SYSTEM_PARAMETERS)
        if not mass_center[5, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 5,
                                   SYSTEM_PARAMETERS)
        if not mass_center[6, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 6,
                                   SYSTEM_PARAMETERS)
        if not mass_center[7, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 7,
                                   SYSTEM_PARAMETERS)
        if not mass_center[8, 6] == 0:
            accel += tc.begin_tree(particles, mass_center, 8,
                                   SYSTEM_PARAMETERS)
    return accel


def tree_code_one_step(particles):
    # Метод leapfrog
    # рассчет изменения компонент скорости, где
    # particles[:, 3] - V_x, Y[:, 4] - V_y, particles[:, 5] - V_z,
    # particles[:, 7] - a_x, particles[:, 8] - a_y, particles[:, 9] - a_z;
    # по формуле V_{i + 1/2} = V_{i} + a_{i} * dt / 2
    particles[:, 3:6] += particles[:, 7:10] * TIME_STEP / 2
    # рассчет изменения координат, где
    # particles[:, 0] - x, particles[:, 1] - y, particles[:, 2] - z
    # по формуле X_{i + 1} = X_{i} + V_{i + 1/2} * dt
    particles[:, 0:3] += particles[:, 3:6] * TIME_STEP
    # Сортируем частицы в порядке возрастания номера ячейки
    particles = cts.distribution(particles, NUMBER_OF_CELLS, CELL_LENGTH)
    # Создаем древоподобную структуру "снизу вверх"
    Cells = cts.create_tree(particles, NUMBER_OF_CELLS, CELL_LENGTH)
    # Проходим древоподобную структуру "сверху вниз"
    # для получения ускорения и потенциала каждой частицы
    # particles[:, 10] - потенциал(phi) на i шаге
    particles[:, 7:11] = tree_root(particles, Cells)
    # Смотрим, есть ли материальные точки, вылетевшие
    # за рассматриваемый объем
    if particles[0, 11] < 0:
        # Считаем взаимодействие с улетевшими точками напрямую
        particles = tc.N_body_direct(particles, SYSTEM_PARAMETERS)
    particles[:, 7:11] *= G
    # рассчет изменения компонент скорости
    # по формуле V_{i + 1} = V_{i + 1/2} + a_{i + 1} * dt / 2
    particles[:, 3:6] += particles[:, 7:10] * TIME_STEP / 2
    return particles


if __name__ == "__main__":
    # v Константы v
    # =======================================================================
    # Гравитационная постоянная
    # G = 6.67408313 * m.pow(10, -11)  # м^3/(кг*с^2)
    G = 4.51811511 * m.pow(10, -7)  # Мпк^3/(М_(Млечного пути)* (10^15 с)^2)
#    G = 1  # Безразмерная константа
    # Скорость света
    # C = 299792458 # м/с
    C = 9.7156188999  # Мпк/(10^15 с)
#    C = 1022.06028776  # Безразмерная константа
# ===========================================================================
# ^ Константы ^
# v Параметры системы v
# ===========================================================================
# Прочие переменные (желательно не трогать)
    MARKER_SIZE = 0.2  # 1
    C_2 = C * C
    ERROR = ''
    ALLOW_LAUNCH = True
    WORKERS = os.cpu_count()

# Временной интервал
    # TIME_STEP = pow(10, 13)  # с
    TIME_STEP = 0.1  # 10^15 с
#    TIME_STEP = 0.1  # Безразмерная величина

# Параметр "сглаживания" гравитационного взаимодействия на близких дистанциях
    SMOOTH_LENGTH = 0.007  # Мпк
#    SMOOTH_LENGTH = 1.0  # Безразмерная величина

# Параметры, которые нужны чаще всего (можно и нужно трогать)
# Количество ячеек по одной оси координат (для tree codes)
# в виде 2^(NUMBER_OF_CELLS)
    NUMBER_OF_CELLS = 2

# Минимальный размер ячейки по одной оси координат
    CELL_LENGTH = SMOOTH_LENGTH * 2000

# Задаем первоначальный размер системы в единицах "CELL_LENGTH"
# для функции parameters_test
    I_CUBE = 3
    J_CUBE = 1
    K_CUBE = 1
    INDENT_I = 0.0
    INDENT_J = 0.0
    INDENT_K = 0.0
    PERIOD = 1.0

# Параметры генерации эллипсоида в единицах (NUMBER_OF_CELLS * CELL_LENGTH / 2)
    A_ELLIPSOID = 1.0
    B_ELLIPSOID = 1.0
    C_ELLIPSOID = 1.0
# Начальные угловые скорости эллипсоида
    W_X = 0.0
    W_Y = 0.0
    W_Z = 0.00154

# Средняя масса наблюдаемых объектов и их пекулярная скорость
    # M_AVG = 1.98892 * pow(10, 41) # кг
    M_AVG = 1  # масс Млечного пути
    # V_AVG = 0 # 4 * pow(10, 5) / np.sqrt(3) # м/с
    V_AVG = 0  # 1.3 * pow(10, -2) / np.sqrt(3) # Мпк/(10^15 c)

# Константы модифицированных теорий гравитации
    ALPHA = 4
    LAMBDA = 3.2407729 * m.pow(10, -29)
    GAMMA = 2.5 * m.pow(10, -6)

# Количество частиц
    N = 10000
# Число шагов
    STEPS = 5000
# Номера шагов, на которых требуется "сфотографировать" положение всех
# материальных точек
    MAKE_PRELAUNCH_SCREENSHOT = True
    SCR_STEP = []
# Тип сгенерированной системы (обязательно заполнить!)
    SYSTEM_GENERATION_TYPE = 'final'
# Тип используемой гравитации
    GRAVITY_TYPE = 'Newton'
# Использовать несколько процессов для вычислений
    USE_MULTIPROCESSING = False
# Использовать данные, введенные вручную
    USE_MANUAL_INPUT = False
# Использовать телеметрию
    USE_TELEMETRY = True
# Обратить время вспять
    INVERSE_TIME = False
# ===========================================================================
# ^ Параметры системы ^
# v Область с исполняемым кодом v
# ===========================================================================
    if USE_MANUAL_INPUT:
        # Осуществляем ввод всех используемых параметров вручную
        INP_PARAMETERS = parameters.manual_input(C)
        print('Делать скриншоты системы?')
        print('y/n')
        INPUT_VARIABLE = input()
        if (INPUT_VARIABLE == 'y') or (INPUT_VARIABLE == 'n'):
            if INPUT_VARIABLE == 'y':
                MAKE_PRELAUNCH_SCREENSHOT = True
                ENABLE_SCREENSHOTS = True
                print('Наберите номера шагов, на которых нужно')
                print('сделать снимок системы')
                print('После того, как все нужные номера введены,')
                print('наберите "end" без кавычек, чтобы продолжить')
                while ENABLE_SCREENSHOTS:
                    INPUT_VAR = input()
                    if INPUT_VAR == 'end':
                        ENABLE_SCREENSHOTS = False
                    else:
                        try:
                            TEMP_VAR = int(INPUT_VAR)
                            SCR_STEP.append(TEMP_VAR)
                        except ValueError:
                            print('Номер шага может быть только целым числом')
            else:
                MAKE_PRELAUNCH_SCREENSHOT = False
                SCR_STEP = []
        else:
            print('Введено недопустимое значение')
            print('Создание скриншотов отменено')
            MAKE_PRELAUNCH_SCREENSHOT = False
        print('Использовать телеметрию?')
        print('y/n')
        INPUT_VARIABLE = input()
        if (INPUT_VARIABLE == 'y') or (INPUT_VARIABLE == 'n'):
            USE_TELEMETRY = INPUT_VARIABLE == 'y'
        else:
            print('Введено недопустимое значение')
            print('Телеметрия не используется')
            USE_TELEMETRY = False
        print('Использовать многоядерность?')
        print('y/n')
        INPUT_VARIABLE = input()
        if (INPUT_VARIABLE == 'y') or (INPUT_VARIABLE == 'n'):
            USE_MULTIPROCESSING = INPUT_VARIABLE == 'y'
        else:
            print('Введено недопустимое значение')
            print('Многоядерность не используется')
        print('Изменить знак у временного интервала?')
        print('y/n')
        INPUT_VARIABLE = input()
        if (INPUT_VARIABLE == 'y') or (INPUT_VARIABLE == 'n'):
            INVERSE_TIME = INPUT_VARIABLE == 'y'
        else:
            print('Введено недопустимое значение')
            INVERSE_TIME = False
#        Достаем из массива значения, необходимые для работы программы
        NUMBER_OF_CELLS = int(INP_PARAMETERS[3])
        TIME_STEP = float(INP_PARAMETERS[18])
        STEPS = int(INP_PARAMETERS[19])
        SMOOTH_LENGTH = float(INP_PARAMETERS[20])
        SYSTEM_GENERATION_TYPE = str(INP_PARAMETERS[21])
        GRAVITY_TYPE = str(INP_PARAMETERS[22])
        ALPHA = float(INP_PARAMETERS[23])
        LAMBDA = float(INP_PARAMETERS[24])
        GAMMA = float(INP_PARAMETERS[25])
    else:
        if NUMBER_OF_CELLS > 0:
            NUMBER_OF_CELLS = int(m.pow(2, int(NUMBER_OF_CELLS)))
        else:
            ALLOW_LAUNCH = False
            print('Количество ячеек не может быть нулевым или отрицательным')
        INP_PARAMETERS = [N, M_AVG, V_AVG,
                          NUMBER_OF_CELLS, CELL_LENGTH,
                          A_ELLIPSOID, B_ELLIPSOID, C_ELLIPSOID,
                          W_X, W_Y, W_Z,
                          I_CUBE, J_CUBE, K_CUBE,
                          INDENT_I, INDENT_J, INDENT_K,
                          PERIOD]
    if (TIME_STEP <= 0) or (CELL_LENGTH <= 0):
        ALLOW_LAUNCH = False
        print('Недопустимые параметры системы')
    if INVERSE_TIME:
        TIME_STEP *= -1
    if (GRAVITY_TYPE == 'plusYukawa') and (LAMBDA == 0):
        ALLOW_LAUNCH = False
        print('Параметр lambda small не может быть нулевым')
    GAMMA = GAMMA * G / C_2
    SYSTEM_PARAMETERS = [NUMBER_OF_CELLS, SMOOTH_LENGTH,
                         GRAVITY_TYPE, ALPHA, LAMBDA, GAMMA]
    try:
        try:
            try:
                # Генерируем или загружаем конфигурацию системы
                # в зваисимости от указанных/введенный параметров
                if (SYSTEM_GENERATION_TYPE == 'cube')\
                    or (SYSTEM_GENERATION_TYPE == 'random')\
                        or (SYSTEM_GENERATION_TYPE == 'ellipsoid'):
                    X = generate.generate_system(SYSTEM_GENERATION_TYPE,
                                                 INP_PARAMETERS)
                    if isinstance(X, str):
                        ALLOW_LAUNCH = False
                        print('Полуоси эллипсоида не могут быть нулевыми')
                    else:
                        np.savetxt('last config.txt', X)
                elif SYSTEM_GENERATION_TYPE == 'last':
                    X = np.loadtxt('last config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif SYSTEM_GENERATION_TYPE == 'debug':
                    X = np.loadtxt('error config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif SYSTEM_GENERATION_TYPE == 'test':
                    X = np.loadtxt('test config.txt', dtype='float64')
                    N = np.size(X, 0)
                elif SYSTEM_GENERATION_TYPE == 'final':
                    X = np.loadtxt('final config.txt', dtype='float64')
                    N = np.size(X, 0)
                else:
                    ALLOW_LAUNCH = False
                    print('Выбранная конфигурация не может быть загружена')
            except IOError:
                ALLOW_LAUNCH = False
                print('Отсутствует необходимый файл конфигурации')
        except TypeError:
            ALLOW_LAUNCH = False
            print('Число материальных точек всегда должно быть целым')
    except ValueError:
        ALLOW_LAUNCH = False
        print('Неприемлимое число материальных точек')
    if ALLOW_LAUNCH:
        # Определяем количество используемых процессов в зависимости
        # от используемого процессора
        if WORKERS >= 8:
            WORKERS = 8
        elif WORKERS >= 4:
            WORKERS = 4
        elif WORKERS >= 2:
            WORKERS = 2
        else:
            USE_MULTIPROCESSING = False
        try:
            if MAKE_PRELAUNCH_SCREENSHOT:
                graph.screenshot_all_system(X, 'Шаг 0', MARKER_SIZE)
            ENERGY = np.zeros([STEPS, 14])
            START_COUNT = time.time()
            for q in range(STEPS):
                # Проверяем, не появились ли в системе сильные
                # гравитационные поля или скорости, равные или
                # превышающие световую
                ERROR = tlm.speed_limit(X, C_2)
                ERROR = tlm.is_gravity_field_weak(X, C_2)
                if not ERROR == '':
                    np.savetxt('error config.txt', X)
                    graph.screenshot_all_system(X, ERROR, MARKER_SIZE)
                    print(ERROR + ' at step ' + str(q))
                    break
                # Выполняем один шаг по алгоритму tree code
                X = tree_code_one_step(X)
                # записываем кинетическую и потенциальную энергию
                # каждой частицы для последующего использования
                ENERGY[q] = tlm.enegry_parameters(q, X)
                X[:, 12] = EM.kinetic_energy_Newton(X)
                X[:, 13] = EM.potential_energy_Newton(X)
                if q == 10000:
                    np.savetxt('backup 10k.txt', X)
                    np.savetxt('backup E 10k.txt', ENERGY)
                if q == 20000:
                    np.savetxt('backup 20k.txt', X)
                    np.savetxt('backup E 20k.txt', ENERGY)
                if q == 30000:
                    np.savetxt('backup 30k.txt', X)
                    np.savetxt('backup E 30k.txt', ENERGY)
                if q in SCR_STEP:
                    graph.screenshot_selected_volume(X, 'Шаг ' + str(q),
                                                     MARKER_SIZE)
            COMPUTING_TIME = time.time() - START_COUNT
            print("Время выполнения", COMPUTING_TIME, "с")
            if USE_TELEMETRY:
                graph.full_telemetry(ENERGY, N)
                graph.plot_combined_energy_in_volume(ENERGY)
                graph.plot_virial_coeff(ENERGY)
        except KeyboardInterrupt:
            print('Работа программы прервана')
            if USE_TELEMETRY:
                graph.full_telemetry(ENERGY, N)
                graph.plot_combined_energy_in_volume(ENERGY)
                graph.plot_virial_coeff(ENERGY)
        print('Сохранить финальную конфигурацию системы?')
        print('y?')
        INPUT_VARIABLE = input()
        if INPUT_VARIABLE == 'y':
            np.savetxt('final config.txt', X)
        else:
            print('Конфигурация не будет сохранена')
        print('Построить графики зависимости M(R)?')
        print('y?')
        if INPUT_VARIABLE == 'y':
            fractal.make_dependency(X, NUMBER_OF_CELLS,
                                    CELL_LENGTH, SMOOTH_LENGTH)
        else:
            print('Построение графиков отменено')
# ===========================================================================
# ^ Область с исполняемым кодом ^
