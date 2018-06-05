# -*- coding: utf-8 -*-

import math


def input_int_value(msg_0, msg_1, msg_2):
    # Функция, обеспечивающая ввод целого положительного числа
    print(msg_0)
    continue_input = True
    while continue_input:
        try:
            variable = int(input())
            if variable > 0:
                continue_input = False
        except ValueError:
            print(msg_1)
            print(msg_2)
    return variable


def input_float_value(msg_0, msg_00, msg_000, msg_1, msg_2):
    # Функция, обеспечивающая ввод рационального положительного числа
    print(msg_0)
    print(msg_00)
    print(msg_000)
    continue_input = True
    while continue_input:
        try:
            variable = float(input())
            if variable >= 0:
                continue_input = False
            else:
                print('Введено некорректное значение. Попробуйте еще раз')
        except ValueError:
            print(msg_1)
            print(msg_2)
    return variable


def input_float_less_1_value(msg_0, msg_00, msg_1, crit):
    # Функция, обеспечивающая ввод рационального числа,
    # не превышающего по модулю 1
    print(msg_0)
    print(msg_00 + str(crit))
    continue_input = True
    while continue_input:
        try:
            variable = float(input())
            if (variable >= -1) and (variable <= 1):
                continue_input = False
            else:
                print('Введено некорректное значение. Попробуйте еще раз')
        except ValueError:
            print(msg_1)
    return variable


def manual_input(c):
    msg_N_0 = 'Введите число материальных точек'
    msg_N_1 = 'Число материальных точек всегда должно быть целым'
    msg_N_2 = 'Введите число материальных точек еще раз'
    msg_n_0 = 'Введите количество ячеек в формате 2^n (нужно задать n)'
    msg_n_1 = 'Число ячеек всегда должно быть целым'
    msg_n_2 = 'Введите число ячеек еще раз'
    msg_steps_0 = 'Введите число временных шагов'
    msg_steps_1 = 'Введено недопустимое число шагов'
    msg_steps_2 = 'Введите число шагов еще раз'
    msg_m_0 = 'Введите среднюю массу материальных точкек в массах галактик'
    msg_m_00 = '(Масса галактики имеет порядок 10^41 кг)'
    msg_m_1 = 'Cредняя масса материальной точки должна быть числом'
    msg_m_2 = 'Введите среднюю массу еще раз'
    msg_v_0 = 'Введите среднюю скорость материальных точкек в Мпк/(10^15 с)'
    msg_v_00 = '(1 Мпк/(10^15 с) = 3,08567758*10^7 м/с)'
    msg_v_000 = 'ВАЖНО ПОМНИТЬ! c = 9.7156188999 Мпк/(10^15 с)'
    msg_v_1 = 'Cредняя скорость материальной точки должна быть числом'
    msg_v_2 = 'Введите среднюю скорость материальных точек еще раз'
    msg_d_0 = 'Введите размер ячейки в Мпк'
    msg_d_1 = 'Размер ячейки должен быть в виде числа'
    msg_d_2 = 'Введите размер ячейки еще раз'
    msg_t_0 = 'Введите временной шаг в единицах (10^15 с)'
    msg_t_1 = 'Временной шаг должен быть в виде числа'
    msg_t_2 = 'Введите временной шаг еще раз'
    msg_ind_0 = 'Введите отступ от границы рассматриваемой'
    msg_ind_i_0 = 'области по оси X в Мпк'
    msg_ind_j_0 = 'области по оси Y в Мпк'
    msg_ind_k_0 = 'области по оси Z в Мпк'
    msg_ind_1 = 'Отступ должен быть в виде числа'
    msg_ind_2 = 'Введите отступ еще раз'
    msg_i_0 = 'Введите число материальных точек по оси X'
    msg_j_0 = 'Введите число материальных точек по оси Y'
    msg_k_0 = 'Введите число материальных точек по оси Z'
    msg_axis_1 = 'Число материальных точек всегда должно быть целым'
    msg_axis_2 = 'Введите число материальных точек еще раз'
    msg_per_0 = 'Введите расстояние между двумя соседними точками,'
    msg_per_00 = 'расположенных на одной оси в единицах длины ячейки'
    msg_per_1 = 'Расстояние должно быть в виде числа'
    msg_per_2 = 'Введите расстояние ячейки еще раз'
    msg_a_0 = 'Введите величину полуоси эллипсоида по оси X'
    msg_b_0 = 'Введите величину полуоси эллипсоида по оси Y'
    msg_c_0 = 'Введите величину полуоси эллипсоида по оси Z'
    msg_abc_0 = 'от 0 до 1. Где 1 соответствует четверти размера системы'
    msg_abc_1 = 'Длина полуоси должна быть числом'
    msg_w_0 = 'Введите начальную угловую скорость в размерности рад/(10^15 с)'
    msg_wx_0 = 'в плоскости YZ. Величина не должна превышать '
    msg_wy_0 = 'в плоскости XZ. Величина не должна превышать '
    msg_wz_0 = 'в плоскости XY. Величина не должна превышать '
    msg_w_1 = 'Угловая скорость должна быть числом'
    msg_eps_0 = 'Введите смягчающую длину потенциала в Мпк'
    msg_eps_1 = 'Смягчающая длина должна быть числом'
    msg_alpha_0 = 'Введите параметр alpha для Ньютоновского потенциала с'
    msg_alpha_1 = 'Параметр alpha должен быть числом'
    msg_alpha_2 = 'Введите параметр alpha еще раз'
    msg_lamb_0 = 'Введите параметр малая lambda для Ньютоновского потенциала с'
    msg_pY_00 = 'добавкой, соответствующей потенциалу Юкавы'
    msg_lamb_1 = 'Параметр малая lambda должен быть числом'
    msg_lamb_2 = 'Введите параметр малая lambda еще раз'
    msg_gamma_0 = 'Введите параметр gamma для потенциала Бранса-Дикке'
    msg_gamma_00 = 'в Ньютоновском пределе'
    msg_gamma_1 = 'Параметр gamma должен быть числом'
    msg_gamma_2 = 'Введите параметр gamma еще раз'
    print('Введите название используемой конфигурации системы')
    # Задаем значения, чтобы программа не выдавала ошибку
    # при загрузке конфигурации
    N = 0
    m_avg, v_avg = 0, 0
    a_inp, b_inp, c_inp = 0, 0, 0
    w_x, w_y, w_z = 0, 0, 0
    i_test, j_test, k_test = 0, 0, 0
    indent_i, indent_j, indent_k = 0, 0, 0
    period = 0
    system_generation_type = str(input())
    Distance = input_float_value(msg_d_0, '', '', msg_d_1, msg_d_2)
    n = input_int_value(msg_n_0, msg_n_1, msg_n_2)
    n = int(math.pow(2, n))
    time_step = input_float_value(msg_t_0, '', '', msg_t_1, msg_t_2)
    Steps = input_int_value(msg_steps_0, msg_steps_1, msg_steps_2)
    eps_smooth = input_float_value(msg_eps_0, '', '', msg_eps_1, '')
    print('Укажите используемую теорию гравитации')
    GRAVITY_TYPE = str(input())
    if not (GRAVITY_TYPE == 'BD') or \
            (GRAVITY_TYPE == 'plusYukawa'):
        GRAVITY_TYPE = 'Newton'
    if GRAVITY_TYPE == 'plusYukawa':
        Alpha = input_float_value(msg_alpha_0, msg_pY_00, '',
                                  msg_alpha_1, msg_alpha_2)
        Lambda = input_float_value(msg_lamb_0, msg_pY_00, '',
                                   msg_lamb_1, msg_lamb_2)
    elif GRAVITY_TYPE == 'BD':
        Gamma = input_float_value(msg_gamma_0, msg_gamma_00, '',
                                  msg_gamma_1, msg_gamma_2)
    else:
        Alpha = 0
        Lambda = 0
        Gamma = 0

    if (system_generation_type == 'random') or \
        (system_generation_type == 'cube') or\
            (system_generation_type == 'ellipsoid'):
        m_avg = input_float_value(msg_m_0, msg_m_00, '', msg_m_1, msg_m_2)
        v_avg = input_float_value(msg_v_0, msg_v_00, msg_v_000,
                                  msg_v_1, msg_v_2)
        if system_generation_type == 'random':
            N = input_int_value(msg_N_0, msg_N_1, msg_N_2)
        elif system_generation_type == 'cube':
            # Ввод отступа "решетки" от границы рассматриваемого объема
            indent_i = input_float_value(msg_ind_0, msg_ind_i_0,
                                         '', msg_ind_1, msg_ind_2)
            indent_j = input_float_value(msg_ind_0, msg_ind_j_0,
                                         '', msg_ind_1, msg_ind_2)
            indent_k = input_float_value(msg_ind_0, msg_ind_k_0,
                                         '', msg_ind_1, msg_ind_2)
            # Ввод количества частиц, приходящегося на каждую ось
            i_test = input_int_value(msg_i_0, msg_axis_1,
                                     msg_axis_2)
            j_test = input_int_value(msg_j_0, msg_axis_1,
                                     msg_axis_2)
            k_test = input_int_value(msg_k_0, msg_axis_1,
                                     msg_axis_2)
            # Ввод расстояния между соседними частицами
            period = input_float_value(msg_per_0, msg_per_00, '',
                                       msg_per_1, msg_per_2)
        elif system_generation_type == 'ellipsoid':
            N = input_int_value(msg_N_0, msg_N_1, msg_N_2)
            w_crit = 2 * c / (n * Distance)
            # Задаем полуоси эллипсоида
            # a_inp - ось x
            # b_inp - ось y
            # c_inp - ось z
            a_inp = input_float_less_1_value(msg_a_0, msg_abc_0,
                                             msg_abc_1, '')
            b_inp = input_float_less_1_value(msg_b_0, msg_abc_0,
                                             msg_abc_1, '')
            c_inp = input_float_less_1_value(msg_c_0, msg_abc_0,
                                             msg_abc_1, '')
            # Задаем начальные угловые скорости материальных точек
            w_x = input_float_less_1_value(msg_w_0, msg_wx_0,
                                           msg_w_1, w_crit)
            w_y = input_float_less_1_value(msg_w_0, msg_wy_0,
                                           msg_w_1, w_crit)
            w_z = input_float_less_1_value(msg_w_0, msg_wz_0,
                                           msg_w_1, w_crit)
    inp_parameters = [N, m_avg, v_avg,
                      n, Distance,
                      a_inp, b_inp, c_inp,
                      w_x, w_y, w_z,
                      i_test, j_test, k_test,
                      indent_i, indent_j, indent_k,
                      period,
                      time_step, Steps, eps_smooth,
                      system_generation_type,
                      GRAVITY_TYPE, Alpha, Lambda, Gamma]
    return inp_parameters
