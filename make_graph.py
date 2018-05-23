# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_max_dE_kinetic(dE):
    # Функция, создающая график максимальной разницы
    # кинетической энергии частиц за все время работы программы
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dE[1:, 0], dE[1:, 4])
    ax.set_xlabel('Номер шага')
    ax.set_ylabel('Kinetic energy')
    ax.set_title('Max kinetic energy difference per step')
    plt.savefig('Максимальное изменение кинетической энергии за шаг', dpi=640)


def plot_max_dE_potential(dE):
    # Функция, создающая график максимальной разницы
    # потенциальной энергии частиц за все время работы программы
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dE[1:, 0], dE[1:, 5])
    ax.set_xlabel('Номер шага')
    ax.set_ylabel('Potential energy')
    ax.set_title('Max potential energy difference per step')
    plt.savefig('Максимальное изменение потенциальной энергии за шаг', dpi=640)


def plot_avg(E, N):
    # Функция, создающая график кинетической энергии частиц
    # за все время работы программы
    Energy = np.copy(E[:, 1:3])
    Energy /= N
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    ax.plot(E[:, 0], Energy[:, 0])
    ax1.plot(E[:, 0], Energy[:, 1])
    ax.xaxis.set_ticklabels([])
    ax1.set_xlabel('Номер шага')
    ax.set_ylabel('Kinetic enegry')
    ax1.set_ylabel('Potential energy')
    ax.set_title('Average energy')
    ax1.set_title(' ')
    plt.savefig('Средняя энергия материальной точки', dpi=640)


def plot_system_enegry(E):
    # Функция, создающая график потенциальной энергии частиц
    # за все время работы программы
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    ax.plot(E[:, 0], E[:, 1])
    ax1.plot(E[:, 0], E[:, 2])
    ax.xaxis.set_ticklabels([])
    ax1.set_xlabel('Номер шага')
    ax.set_ylabel('Kinetic enegry')
    ax1.set_ylabel('Potential energy')
    ax.set_title('Energy at step')
    plt.savefig('Кинетическая и потенциальная энергия системы', dpi=640)


def plot_total_energy(E):
    # Функция, создающая график потенциальной энергии частиц
    # за все время работы программы
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(E[:, 0], E[:, 3])
    ax.set_xlabel('Номер шага')
    ax.set_ylabel('Энергия')
    ax.set_title('Полная энергия системы')
    plt.savefig('Полная энергия системы', dpi=640)


def plot_combined_energy(E):
    # Функция, создающая график потенциальной энергии частиц
    # за все время работы программы
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(E[:, 0], E[:, 3], label='Полная энергия', color='black')
    ax.plot(E[:, 0], E[:, 1], label='Кинетическая энергия', color='red')
    ax.plot(E[:, 0], E[:, 2], label='Потенциальная энергия', color='blue')
    ax.set_xlabel('Номер шага')
    ax.set_ylabel('Энергия')
    ax.set_title('Полная энергия системы')
    plt.legend()
    plt.savefig('Кинетическая, потенциальная, полная энергия системы', dpi=640)


def screenshot(System_parameters, name, point_size):
    # Функция для "скирншота" положения всех частиц
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = System_parameters[:, 0]
    y = System_parameters[:, 1]
    z = System_parameters[:, 2]
    ax.scatter(x, y, z, color='red', s=point_size)
    ax.autoscale(False)
    ax.set_xlabel('x, кпк')
    ax.set_ylabel('y, кпк')
    ax.set_zlabel('z, кпк')
    plt.savefig(name, dpi=1280)
#    plt.show()


def full_telemetry(E, N):
    plot_max_dE_kinetic(E)
    plot_max_dE_potential(E)
    plot_avg(E, N)
    plot_system_enegry(E)
    plot_total_energy(E)
    plot_combined_energy(E)
