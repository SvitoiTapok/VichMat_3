import tkinter as tk
from tkinter import ttk,filedialog, messagebox

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import sympy as sp
from numpy.lib.format import read_magic
from sympy.physics.units import action


def left_pram_meth(func, left, right, n, acc):
    x = []
    y = []
    h = (right-left)/n
    ans = 0
    for i in range(n):
        ar = left+i*h
        zn = careful_solve(func, ar, acc, left, right, h)
        ans += zn
        x.append(ar)
        y.append(zn/h)
    return [ans, x, y]

def right_pram_meth(func, left, right, n, acc):
    x = []
    y = []
    h = (right-left)/n
    ans = 0
    for i in range(1, n+1):
        ar = left + i * h
        zn = careful_solve(func, ar, acc, left, right, h)
        ans += zn
        x.append(ar)
        y.append(zn/h)
    return [ans, x, y]

def mid_pram_meth(func, left, right, n, acc):
    x = []
    y = []
    h = (right-left)/n
    ans = 0
    for i in range(n):
        ar = left + (i+1/2) * h
        zn = careful_solve(func, ar, acc, left, right, h)
        ans += zn
        x.append(ar)
        y.append(zn/h)
    return [ans, x, y]

def trap_meth(func, left, right, n, acc):

    h = (right-left)/n
    ans = (careful_solve(func, left, acc, left, right, h) + careful_solve(func, right, acc, left, right, h))/2*h
    x = [left, right]
    y = [careful_solve(func, left, acc, left, right, h), careful_solve(func, right, acc, left, right, h)]
    for i in range(1, n):
        ar = left + i * h
        zn = careful_solve(func, ar, acc, left, right, h)
        ans += zn
        x.append(ar)
        y.append(zn/h)
    return [ans, x, y]

def simps_meth(func, left, right, n, acc):

    h = (right-left)/n
    ans = (careful_solve(func, left, acc, left, right, h) + careful_solve(func, right, acc, left, right, h))/3*h
    x = [left, right]
    y = [careful_solve(func, left, acc, left, right, h), careful_solve(func, right, acc, left, right, h)]
    for i in range(1, n):
        ar = left + i * h
        zn = careful_solve(func, ar, acc, left, right, h)
        x.append(ar)
        y.append(zn/h)
        if(i%2):
            ans += 4*zn/3
        else:
            ans += 2*zn/3
    return [ans, x, y]

def careful_solve(func, ar, acc, left, right, h):
    acc /= 100
    try:
        zn = func(ar) * h
    except:
        if ar == left:
            zn = func(ar + acc) * h
        elif ar == right:
            zn = func(ar + acc) * h
        else:
            zn = (func(ar - acc) + func(ar + acc)) / 2 * h
    return zn


def func1(x):
    return np.sin(np.exp(x) + x)


def func2(x):
    return x ** np.sin(x) - 0.5 * x


def func3(x):
    return np.arctan(x) * np.sin(x) * 9


def func4(x):
    return x ** 3 - 9 * x ** 2 + x + 11


def func5(x):
    return -3*x**3-5*x**2+4*x-2

def func6(x):
    return 1/abs(x)**0.5

def func7(x):
    return 1/(1-x)

# print(left_pram_meth(func5, -3, -1, 1000))
# print(right_pram_meth(func5, -3, -1, 1000))
# print(mid_pram_meth(func5, -3, -1, 10))
# print(trap_meth(func5, -3, -1, 10))
# print(simps_meth(func5, -3, -1, 10))

def draw_func():
    func = cur_func
    name = cur_func_name
    ax.clear()
    try:
        left, right = float(left_gran.get().replace(',', '.')), float(right_gran.get().replace(',', '.'))
        if left>=right:
            messagebox.showerror("Ошибка!", "Сообщение об ошибке: левая граница больше правой")
            return 1
        x = np.linspace(left, right, 100)
        y = func(x)
    except:
        messagebox.showerror("Ошибка!", "Сообщение об ошибке: некорректные значения границ")
        return 1
    ax.plot(x, y)
    ax.set_title(name)
    ax.grid( color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.axhline(0, color='black', linewidth=0.5)
    canvas.draw()
    return 0

def calculate_int():
    delete_pre_res_func()
    if draw_func():
        return 1
    meth = cur_meth
    name = cur_meth_name
    func = cur_func
    left, right = float(left_gran.get().replace(',', '.')), float(right_gran.get().replace(',', '.'))
    try:
        acc = float(accuracy.get().replace(',', '.'))
        if acc<=0: raise Exception
    except:
        messagebox.showerror("Ошибка!", "Сообщение об ошибке: некорректное значение точности")
        return 1
    #calc_err_label.config(text="Сообщение об ошибке:")
    n=8
    I1 = meth(func, left, right, n//2, acc)[0]
    I2 = meth(func, left, right, n, acc)[0]
    k=2
    if name=="Метод Симпсона":
        k=4
    while abs(I2-I1)/(2**k-1)>acc:

        n*=2
        I1 = I2
        I2 = meth(func, left, right, n, acc)[0]
        print(n)
        if n>10000000:
            messagebox.showerror("Ошибка!",f"Сообщение об ошибке: Интеграл расходящийся(за {n} разбиений сходимости не обнаружено)")
            return 1

    answ = meth(func, left, right, n, acc)
    #ax.scatter(answ[1], answ[2], color='red', label='точки рассчета функции', s=30)
    canvas.draw()
    Info_label.config(text=f"Информация о выполнении({name})")
    root_label.config(text=f"Рассчитанное значение интеграла: {round(I2, 10)}")
    iter_label.config(text=f"Число разбиений: {n}")



def on_func_combobox_change(event):
    delete_pre_res_func()
    # Функция, которая вызывается при изменении значения в выпадающем списке
    global cur_func
    global cur_func_name
    selected_value = func_combobox.get()  # Получаем выбранное значение
    if selected_value=="sin(exp(x) + x)":
        cur_func=func1
    elif selected_value=="x ** sin(x) - 0.5 * x":
        cur_func=func2
    elif selected_value=="atan(x) * sin(x) * 9":
        cur_func=func3
    elif selected_value=="x ** 3 - 9 * x ** 2 + x + 11":
        cur_func=func4
    elif selected_value=="-3*x**3 - 5*x**2 + 4*x - 2":
        cur_func=func5
    elif selected_value=="1 / abs(x) ** 0.5":
        cur_func=func6
    elif selected_value=="1 / (1 - x)":
        cur_func=func7
    cur_func_name = selected_value
    func_combobox_label.config(text=f"Вы выбрали: {selected_value}")  # Обновляем текст метки

def on_meth_combobox_change(event):
    delete_pre_res_func()
    global cur_meth
    global cur_meth_name
    selected_value = method_combobox.get()

    if selected_value=="Метод левых прямоугольников":
        cur_meth=left_pram_meth
    elif selected_value=="Метод правых прямоугольников":
        cur_meth=right_pram_meth
    elif selected_value=="Метод средних прямоугольников":
        cur_meth=mid_pram_meth
    elif selected_value=="Метод трапеций":
        cur_meth=trap_meth
    elif selected_value=="Метод Симпсона":
        cur_meth=simps_meth

    cur_meth_name = selected_value
    method_combobox_label.config(text=f"Вы выбрали: {selected_value}")  # Обновляем текст метки

def delete_pre_res_func():
    Info_label.config(text="Информация о выполнении")
    root_label.config(text=f"Рассчитанное значение интеграла:")
    iter_label.config(text=f"Число разбиений:")

# def delete_pre_res_sys():
#     root_label1.config(text=f"Рассчитанный конень:")
#     iter_label1.config(text=f"Количество итераций:")
#     func_label1.config(text=f"Значение в точке:")


def show_frame(frame):
    # Скрыть все фреймы
    # for f in (int_frame, nes_frame):
    #     f.grid_forget()
    # Показать выбранный фрейм
    frame.grid(row=0, column=0, sticky="nsew")

def save():
    if root_label.cget('text')=="Рассчитанное значение интеграла:":
        messagebox.showerror("Ошибка!",
                             f"Сообщение об ошибке: Перед сохранением рассчитайте интеграл, информацию о рассчете которого хотите сохранить")
        return 1
    filepath = filedialog.asksaveasfilename(
        defaultextension=".txt",
        filetypes=[("Текстовые файлы", "*.txt")]
    )

    if not filepath:
        return 1

    with open(filepath, "w", encoding="utf-8") as file:
        file.write(f"Функция:{cur_func_name}, Метод:{cur_meth_name}, Точность:{accuracy.get()}, Левая граница:{left_gran.get()}, Правая граница:{right_gran.get()}, {root_label.cget('text')}, {iter_label.cget('text')}\n")
    messagebox.showinfo("Успех", f"Файл сохранён: {filepath}")

cur_func=func1
cur_func_name="sin(exp(x) + x)"
cur_meth=left_pram_meth
cur_meth_name="Метод левых прямоугольников"

root = tk.Tk()
root.title("Динамическое построение графиков")
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

int_frame = ttk.Frame(root)
int_frame.grid_rowconfigure(0, weight=1)
int_frame.grid_rowconfigure(1, weight=1)
int_frame.grid_rowconfigure(2, weight=1)
int_frame.grid_rowconfigure(3, weight=1)
int_frame.grid_rowconfigure(4, weight=1)
int_frame.grid_rowconfigure(5, weight=1)
int_frame.grid_rowconfigure(6, weight=1)
int_frame.grid_rowconfigure(7, weight=1)
int_frame.grid_rowconfigure(8, weight=1)
int_frame.grid_rowconfigure(9, weight=1)
int_frame.grid_rowconfigure(10, weight=1)
int_frame.grid_columnconfigure(0, weight=1)
int_frame.grid_columnconfigure(1, weight=1)

# nes_frame = ttk.Frame(root)
# nes_frame.grid_rowconfigure(0, weight=1)
# nes_frame.grid_rowconfigure(1, weight=1)
# nes_frame.grid_rowconfigure(2, weight=1)
# nes_frame.grid_rowconfigure(3, weight=1)
# nes_frame.grid_rowconfigure(4, weight=1)
# nes_frame.grid_rowconfigure(5, weight=1)
# nes_frame.grid_rowconfigure(6, weight=1)
# nes_frame.grid_rowconfigure(7, weight=1)
# nes_frame.grid_rowconfigure(8, weight=1)
# nes_frame.grid_rowconfigure(9, weight=1)
# nes_frame.grid_rowconfigure(10, weight=1)
# nes_frame.grid_rowconfigure(11, weight=1)
# nes_frame.grid_rowconfigure(12, weight=1)
# nes_frame.grid_rowconfigure(13, weight=1)
# nes_frame.grid_columnconfigure(0, weight=1)
# nes_frame.grid_columnconfigure(1, weight=1)

# Поле для ввода погрешности
accuracy_label = ttk.Label(int_frame, text="Точность:")
accuracy_label.grid(row=0,column=0, padx=10, pady=10,sticky="ew")
accuracy = ttk.Entry(int_frame)
accuracy.insert(0, "0.01")
accuracy.grid(row=0,column=1, padx=10, pady=10,sticky="ew")

# Поле для ввода отрезка
left_gran_label = ttk.Label(int_frame, text="левая граница:")
left_gran_label.grid(row=1,column=0, padx=10, pady=10,sticky="ew")
left_gran = ttk.Entry(int_frame)
left_gran.insert(0, "0.0")
left_gran.grid(row=1,column=1, padx=10, pady=10,sticky="ew")


right_gran_label = ttk.Label(int_frame, text="правая граница:")
right_gran_label.grid(row=2,column=0, padx=10, pady=10,sticky="ew")
right_gran = ttk.Entry(int_frame)
right_gran.insert(0, "1.0")
right_gran.grid(row=2,column=1, padx=10, pady=10,sticky="ew")
functions = ["sin(exp(x) + x)", "x ** sin(x) - 0.5 * x", "atan(x) * sin(x) * 9", "x ** 3 - 9 * x ** 2 + x + 11", "-3*x**3 - 5*x**2 + 4*x - 2", "1 / abs(x) ** 0.5", "1 / (1 - x)"]
func_combobox = ttk.Combobox(int_frame, values=functions)
func_combobox.set("sin(exp(x) + x)")
#func_combobox.pack(pady=10)
func_combobox.bind("<<ComboboxSelected>>", on_func_combobox_change)
func_combobox.grid(row=3,column=0)

# какая функция выбрана
func_combobox_label = tk.Label(int_frame, text="Выберите вариант функции из списка")
func_combobox_label.grid(row=4,column=0)


fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=int_frame)
draw_button = ttk.Button(int_frame, text="Нарисовать график", command=draw_func)
err_label = ttk.Label(int_frame, text="Сообщение об ошибке:")
draw_button.grid(row=5, column=0, padx=10, pady=10,sticky="ew")
err_label.grid(row=6, column=0, padx=10, pady=10,sticky="ew")

methods = ["Метод левых прямоугольников", "Метод правых прямоугольников", "Метод средних прямоугольников", "Метод трапеций", "Метод Симпсона"]
method_combobox = ttk.Combobox(int_frame, values=methods)
method_combobox.set("Метод левых прямоугольников")
#func_combobox.pack(pady=10)
method_combobox.bind("<<ComboboxSelected>>", on_meth_combobox_change)
method_combobox.grid(row=3,column=1)

# какой метод выбрана
method_combobox_label = tk.Label(int_frame, text="Выберите метод из списка")
method_combobox_label.grid(row=4,column=1)

calc_button = ttk.Button(int_frame, text="Расчитать интеграл", command=calculate_int)
calc_err_label = ttk.Label(int_frame, text="Сообщение об ошибке:")
calc_button.grid(row=5, column=1, padx=10, pady=10,sticky="ew")
calc_err_label.grid(row=6, column=1, padx=10, pady=10,sticky="ew")



save_button = ttk.Button(int_frame, text="Сохранить результат", command=save)
save_button.grid(row=7, column=0, columnspan=2)
canvas.get_tk_widget().grid(row=8, column=0, columnspan=2)

Info_label = ttk.Label(int_frame, text="Информация о выполнении")
Info_label.grid(row=9, column=0, columnspan=2)


root_label = ttk.Label(int_frame, text="Рассчитанное значение интеграла:")
root_label.grid(row=10, column=0, padx=10, pady=10,sticky="ew")

iter_label = ttk.Label(int_frame, text="Число разбиений:")
iter_label.grid(row=10, column=1, padx=10, pady=10,sticky="ew")

#func_combobox_label.pack(pady=10)
################################################################################################################################################

# accuracy_label = ttk.Label(nes_frame, text="Точность:")
# accuracy_label.grid(row=0,column=0, padx=10, pady=10,sticky="ew")
# accuracy1 = ttk.Entry(nes_frame)
# accuracy1.insert(0, "0.01")
# accuracy1.grid(row=0,column=1, padx=10, pady=10,sticky="ew")
#
# # Поля для ввода отрезков
# left_granx_label = ttk.Label(nes_frame, text="левая граница x:")
# left_granx = ttk.Entry(nes_frame)
# left_granx.grid(row=1,column=1, padx=10, pady=10,sticky="ew")
# left_granx.insert(0, "0.0")
# left_granx_label.grid(row=1,column=0, padx=10, pady=10,sticky="ew")
#
#
# right_granx_label = ttk.Label(nes_frame, text="правая граница x:")
# right_granx = ttk.Entry(nes_frame)
# right_granx.insert(0, "1.0")
# right_granx.grid(row=2,column=1, padx=10, pady=10,sticky="ew")
# right_granx_label.grid(row=2,column=0, padx=10, pady=10,sticky="ew")
#
#
# left_grany_label = ttk.Label(nes_frame, text="нижняя граница y:")
# left_grany = ttk.Entry(nes_frame)
# left_grany.grid(row=3,column=1, padx=10, pady=10,sticky="ew")
# left_grany.insert(0, "0.0")
# left_grany_label.grid(row=3,column=0, padx=10, pady=10,sticky="ew")
#
#
# right_grany_label = ttk.Label(nes_frame, text="верхняя граница y:")
# right_grany = ttk.Entry(nes_frame)
# right_grany.insert(0, "1.0")
# right_grany.grid(row=4,column=1, padx=10, pady=10, sticky="ew")
# right_grany_label.grid(row=4,column=0, padx=10, pady=10, sticky="ew")
#
#
#
# systems = ["0.1*x**2+x+0.2*y**2-0.3      0.2*x**2 + y + 0.1*x*y-0.7",
#            "log(y + 2) - x      exp(-x) - y",
#            "(y+np.exp(-x))/3-x      (x**2+1)/4-y"]
# sys_combobox = ttk.Combobox(nes_frame, values=systems)
# sys_combobox.set("0.1*x**2+x+0.2*y**2-0.3\n0.2*x**2 + y + 0.1*x*y-0.7")
# sys_combobox.grid(row=5,column=0, padx=10, pady=10,sticky="ew")
# # sys_combobox.bind("<<ComboboxSelected>>", on_sys_combobox_change)
#
#
# # какая система выбрана
# sys_combobox_label = tk.Label(nes_frame, text="Выберите вариант системы из списка")
# sys_combobox_label.grid(row=6,column=0, padx=10, pady=10,sticky="ew")
#
#
#
#
#
#
# # fig1, ax1 = plt.subplots()
# # canvas1 = FigureCanvasTkAgg(fig1, master=nes_frame)
# # draw_button1 = ttk.Button(nes_frame, text="Нарисовать множество точек", command=draw_sys)
# # err_label1 = ttk.Label(nes_frame, text="Сообщение об ошибке:")
# # draw_button1.grid(row=7, column=0, padx=10, pady=10,sticky="ew")
# # err_label1.grid(row=8, column=0, padx=10, pady=10,sticky="ew")
# # switch_button_func = ttk.Button(nes_frame, text="Перейти к функциям", command=lambda: show_frame(int_frame))
# # switch_button_func.grid(row=9,column=0, columnspan=2)
# # canvas1.get_tk_widget().grid(row=10, column=0, columnspan=2)
# #
# # calc_button1 = ttk.Button(nes_frame, text="Расчитать корень", command=calculate_root_sys)
# # calc_err_label1 = ttk.Label(nes_frame, text="Сообщение об ошибке:")
# # calc_button1.grid(row=6, column=1, padx=10, pady=10,sticky="ew")
# # calc_err_label1.grid(row=7, column=1, padx=10, pady=10,sticky="ew")
# #
# # Info_label1 = ttk.Label(nes_frame, text="Информация о выполнении")
# # Info_label1.grid(row=11, column=0, columnspan=2)
# #
# # root_label1 = ttk.Label(nes_frame, text="Рассчитанный конень:")
# # root_label1.grid(row=12, column=0, padx=10, pady=10,sticky="ew")
# #
# # iter_label1 = ttk.Label(nes_frame, text="Количество итераций:")
# # iter_label1.grid(row=12, column=1, padx=10, pady=10,sticky="ew")
# #
# # func_label1 = ttk.Label(nes_frame, text="Значение в точке:")
# # func_label1.grid(row=13, column=0, columnspan=2)
#
#
#
# # # Кнопка для обновления графика
#
# # update_button.grid(row=2, column=0, columnspan=2, pady=10)
#
# # Создание графика
#
#
#
#
# #canvas.get_tk_widget().pack()
# #field_set_1=[accuracy_label, accuracy, left_gran_label, left_gran, right_gran_label, right_gran, func_combobox, func_combobox_label, op_combobox, draw_button, canvas.get_tk_widget(), err_label]
# #field_set_2=[accuracy_label, accuracy, left_granx_label, left_granx, rigth_granx_label, right_granx, left_grany_label, left_grany, right_grany_label, right_grany, sys_combobox, sys_combobox_label, op_combobox, canvas.get_tk_widget(), err_label]
# #for widget in field_set_1:
# #    widget.pack(pady=5)
#
#
# # Инициализация графика
# # update_plot()
show_frame(int_frame)
# Запуск основного цикла
root.mainloop()
