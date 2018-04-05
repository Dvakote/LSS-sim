# LSS - simulation
Этот readme предназначен для N-body.py от 05.04.2018.
Изменения:
В двух словах: исправлена целая серия багов, которые присутствовали
в программе изначально, но были не заметны.

Подробно:
Теперь взаимодействие частица-частица обрабатывается корректно.
Номера подъячеек теперь хранятся в той же строке, где и 
параметры центра масс.
Исправлен баг, из-за которого не все ячейки обрабатывались.
Изменен порядок исполнения для части кода, определяющего
количество пройденных ячеек.

======================================================

N-body.py от 03.04.2018
Изменения:
В двух словах: переработан алгоритм, отвечающий за распределение
частиц по ячейкам и вычисление параметров центров масс. В результате
существенно повысилось быстродействие всего алгоритма Tree code в
рассчете на один шаг.

Подробно:
Важные изменения:
Полностью переписана функция Tree_code_gravity(Y), разделена работа с
частицами и работа с ячейками, изменен порядок получения параметров ячеек.
Теперь параметры центра масс только самых малых ячеек находятся исходя из
параметров частиц в ячейке (новая функция Particles_to_cell.
Параметры центров масс для остальных ячеек теперь находятся из уже вложенных
ячеек (новая функция Cells_to_cell).

Изменена функция для распределения частиц по ячейкам Distribution,
теперь в она еще и сортирует частицы по номерам ячеек. Сортировка 
осуществляется встроенным в пакет numpy методом (сложность алгоритма O(N logN)).

Небольшие изменения:
В параметры частиц при генерации randomize_parameters() заранее
включают в себя номер ячейки, где частица находится (по-умолчанию "0")

Название функций для получения матрицы с параметрами частиц было изменено
с spawn_test и spawn_random на birth_test и birth_random соответственно.

В матрицу, где находятся параметры ячеек (R_final) на 4 и 5 позиции
перемещенны номера парвой частицы в ячейке и последней, соответственно.
С 7 на 3 позицию перемещена степень двойки, указывающая на число 
ячеек по одной оси.

Переработены функции Calculate_cell и Calculate_cube под новый формат
данных (включение матрицы cell_count в R_final и переход от координат
ячеек к порядковому номеру)

======================================================

N-body.py от 28.03.2018.
Так как я проводил тестирование через изменение переменных
внутри программы, то дата создания отличается от указанной,
но по структуре и функциональности они идентичны.

Сама программа написана на языке Python 3.6.3 в сборке 
Anaconda3 5.0.1, но скорее всего запустится при наличии
всех необходимых пакетов.

Большая часть комментариев касательно работы той или иной 
части уже есть в коде программы, но если есть какие то 
непонятные моменты, то спрашивайте.

Текст программы разделен на секции в виде:

# v "Название секции" v
#===========================================================================
... "Содержание" ...
#===========================================================================
# ^ "Название секции" ^

надеюсь, это поможет упростить понимание кода.
Список секций:
1) Подключаемые пакеты
2) Константы
3) Параметры системы
4) Используемые функции
5) Область с исполняемым кодом

Думаю, комментировать первую секцию будет излишним.

Во второй секции собраны константы, отвечающие за физику,
значения приблизительные, хотя и подбирались как 
приближенные к реальным.
В закомментированном виде те же значения, что и используемые
для вычислений, но уже в системе СИ.

Третья секция содержит константы, отвечающие за работу
программы. Здесь можно указать количество рассматриваемых
частиц, количество шагов, которое должна выполнить система,
а так же начальный объем рассматриваемой системы.
! Размер системы считается как число ячеек (n) * минимальный 
размер ячейки (Distance) !
!! Если использовать Tree code, то количество ячеек (n) будет
напрямую влиять на производительность и не всегда в 
положительную сторону !!
!!! Указанные параметры работают только для случайной генерации,
если нужно проверить конкретный случай то нужно изменять 
функцию parameters_test. !!!

В четвертой секции содержаться все используемые функции.
Обычно распологаются в порядке вызова в основной части 
программы, но если функция вызывает другую функцию, то 
вызываемая обязательно находится "выше". Исключением
является рекурсивная часть.

Пятая секция отведена под исполняемый код. Здесь уже есть
как вычисление сил напрямую, так и алгоритмом Tree code,
нужно только раскомментировать нужную часть и закомментировать
другую. Если нужна информация по общему поведению частиц,
но не хочется разбирать все это в виде матрицы, то можно сделать
"скриншоты" положения всех частиц в системе. Так же можно 
замерить время выполнения программы или ее части, если 
раскомментировать нужные функции.

В данный момент метод Tree code отработан только на 1 шаг,
в принципе, ничего не мешает сделать N шагов, но за стабильность
не ручаюсь. Так же в этом алгоритме есть две функции, которые
мне так и не удалось оттестировать полностью, но результат,
вроде, в пределах нормы.
