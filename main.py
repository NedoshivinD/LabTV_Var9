from tkinter import *
import random
from tkinter import ttk
from PIL import Image, ImageTk
import numpy as np
from collections import Counter
from operator import itemgetter
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib
import pandas as pd

window = Tk()
window.title('Лабораторная работа №1')
window.geometry("1200x780")

arr_sum = []
Arr_F = []
p_num = []
xs = []


def task_window():
    task = Toplevel()
    task.title('Задача')
    img = ImageTk.PhotoImage(Image.open("./task.png"))
    panel = Label(task, image=img)
    panel.pack(side='bottom', fill='both', expand='yes')
    task.mainloop()

def graphs_window():
    graphs = Toplevel()
    graphs.title('Графики')
    f = Figure(figsize=(5, 5), dpi=100)
    x = []
    r = int(r_value.get())
    x.append(0)
    for i in range(3, 3*r+2):
        x.append(i)
    f = Figure(figsize=(5, 5), dpi=100)
    a = f.add_subplot(111)

    a.step(x, Arr_F, "g-o", where='pre', label='Теоретическая')
    a.step(xs, p_num, "b-o", where='pre', label='Выборочная')
    canvas = FigureCanvasTkAgg(f, graphs)
    canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, graphs)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
    a.legend()

    graphs.mainloop()

def run():
    r = int(r_value.get())
    N1 = int(N_value.get())
    ver1 = float(ver1_value.get())
    ver2 = float(ver2_value.get())
    ver3 = float(ver3_value.get())
    K = int(k_value.get())
    alp = float(alpha_value.get())
    arr_sum = []
    x_i = 0

    for i in range(0, N1):
        sum1 = 1
        sum2 = 1
        sum3 = 1
        flag1 = 0
        flag2 = 0
        flag3 = 0
        for j in range(0, r-1):
            p1 = random.random()
            p2 = random.random()
            p3 = random.random()
            if (p1 > ver1 and flag1==0):
                sum1 = sum1 + 1
            else:
                flag1 = 1
            if (p2 > ver2 and flag2==0):
                sum2 = sum2 + 1
            else:
                flag2 = 1
            if (p3 > ver3 and flag3==0):
                sum3 = sum3 + 1
            else:
                flag3 = 1
        sum = sum1 + sum2 + sum3
        arr_sum.insert(i, sum)

    c = Counter(arr_sum)
    c.most_common()
    c = sorted(c.items())
    refresh_table1(c)
    arr_sum.sort()
    Arr_sum = []
    for i in range(r, 3 * r + 1):
        Arr_sum.append(i)
    lst = []
    for i in range(len(c)):
        lst.append(c[i][0])
    result = list(set(Arr_sum) - set(lst))
    ni = []
    j = 0
    if lst[0] != Arr_sum[0]:
        ni.append(0)
        j += 1
    i = 0
    while len(lst) != len(Arr_sum):
        lst.append(0)
    while j < len(Arr_sum) and i < len(lst):
        if (Arr_sum[j] == lst[i]):
            ni.append(c[i][1])
            i += 1
            j += 1
        else:
            ni.append(0)
            j += 1
    list_ = []
    ver1_ = []
    ver2_ = []
    ver3_ = []
    mat_og1 = 0
    mat_og2 = 0
    mat_og3 = 0
    mat_og1_ = 0
    mat_og2_ = 0
    mat_og3_ = 0
    for i in range(1, r + 1):
        ver_1 = ((1 - ver1) ** (i - 1)) * ver1
        ver_2 = ((1 - ver2) ** (i - 1)) * ver2
        ver_3 = ((1 - ver3) ** (i - 1)) * ver3
        if (i == r):
            ver_1 += (1 - ver1) ** r
            ver_2 += (1 - ver2) ** r
            ver_3 += (1 - ver3) ** r
        ver1_.append(ver_1)
        ver2_.append(ver_2)
        ver3_.append(ver_3)
        mat_og1 += ver_1 * i
        mat_og2 += ver_2 * i
        mat_og3 += ver_3 * i
        mat_og1_ += ver_1 * (i ** 2)
        mat_og2_ += ver_2 * (i ** 2)
        mat_og3_ += ver_3 * (i ** 2)
    mat_og = mat_og1 + mat_og2 + mat_og3
    mat_og_ = mat_og1_ + mat_og2_ + mat_og3_
    list_.append(mat_og)
    for i in range(len(arr_sum)):
        x_i += arr_sum[i]
    x_ = x_i / N1

    list_.append(x_)
    abs_ = abs(mat_og - x_)
    list_.append(abs_)
    disp = mat_og_ + 2 * mat_og1 * mat_og2 + 2 * mat_og1 * mat_og3 + 2 * mat_og2 * mat_og3 - mat_og ** 2

    list_.append(disp)
    S = 0
    for i in range(len(arr_sum)):
        S += (arr_sum[i] - x_) ** 2
    S2 = S / N1

    list_.append(S2)

    abs__ = abs(disp - S2)

    list_.append(abs__)
    if (N1 % 2 == 1):
        k = int((N1 - 1) / 2)
        Me = arr_sum[k]
    else:
        k = int(N1 / 2)
        Me = (arr_sum[k - 1] + arr_sum[k]) / 2

    list_.append(Me)

    R = arr_sum[-1] - arr_sum[0]

    list_.append(R)
    refresh_table2(list_)

    vers = []
    for k in range(3, 3 * r + 1):
        vers_ = 0
        for i in range(1, r+1):
            for j in range(1, r+1):
                if ((k - i - j - 1) >= 0 and (k - i - j) <= len(ver3_)):
                    vers_ += ver1_[i-1] * ver2_[j-1] * ver3_[k - i - j - 1]
        vers.append(round(vers_, 9))
    sums = 0
    Arr_F.clear()
    Arr_F.append(0)
    Arr_F.append(0)
    for i in range(len(vers)):
        sums += vers[i]
        Arr_F.append(sums)

    p_num_ = []
    p_num.clear()
    sum = 0
    p_num.insert(0, 0)
    p_num.insert(0, 0)
    for i in range(len(c)):
        sum += float(c[i][1]/N1)
        p_num.append(sum)
    sum = 0
    p_num_.insert(0, 0)
    p_num_.insert(0, 0)
    for i in range(len(c)):
        sum += float(c[i][1]/N1)
        p_num_.append(sum)
    xs.clear()

    znah = []
    for i in range(len(c)):
        znah.append(int(c[i][0]))
    for i in range(len(c)):
        xs.append(int(c[i][0]))
    p_num.sort()
    for i in range(r, 3 * r + 1):
        if i not in znah:
            xs.append(i)
            if i < min(znah):
                p_num.append(0)
            elif i > max(znah):
                p_num.append(1)
            else:
                p_num.append(p_num_[len(p_num_) - i - 1])

    xs.sort()
    xs.append(3 * r + 1)
    xs.insert(0, 0)
    p_num.sort()

    r_ = []
    for i in range(len(p_num)):
        r_.append(abs(p_num[i] - Arr_F[i]))
    D = max(r_)

    D_ = Label(text='Мера расхождения: ' + 'D = ' + str(round(D, 7)))
    D_.grid(row=10, column=0, columnspan=2, padx=10, pady=0, sticky="we")

    refresh_table3(c, vers)

    r__ = []
    for i in range (len(c[1])):
        r__.append(abs(float(c[i][1]/N1) - vers[i]))

    maxx = max(r__)

    R_ = Label(text='max|nj / n - P(yj)| = ' + str(round(maxx, 7)))
    R_.grid(row=11, column=0, columnspan=2, padx=10, pady=0, sticky="we")

    #3 часть
    def Gamma_func(value):
        if (value == 1):
            return value
        if (value == 0.5):
            return np.sqrt(np.pi)
        return ((value - 1) * Gamma_func(value - 1))
    def intervals(k):
        tmp = vers
        num = len(tmp) / k
        t = []
        d = []
        ntmp = []
        nj_ = []
        nj = []
        count = 0
        count_2 = num
        for i in Arr_sum:
            t.append(i)
            ntmp.append(ni[count])
            count += 1
            if (count >= count_2):
                d.append(t)
                nj_.append(ntmp)
                ntmp = []
                t = []
                count_2 += num
        for i in nj_:
            nj.append(np.sum(i))
        return d, nj

    def q(d):
        q = []
        for i in d:
            q.append(np.sum(i))
        return q

    def R0(qj, nj):
        summa = 0
        for i in range(len(qj)):
            summa += (nj[i] - N1 * qj[i]) ** 2 / (N1 * qj[i])
        return summa

    def fx(x, r):
        if (x > 0):
            return (2 ** (-r / 2)) * (Gamma_func(r / 2) ** (-1)) * (x ** (r / 2 - 1) * (np.exp(-x / 2)))
        else:
            return 0

    def F_(R0, k, n=1000):
        a = 0
        b = R0
        summa = 0
        k = k - 1
        for i in range(1, n + 1, 1):
            summa += (b - a) * (fx(a + (b - a) * ((i - 1) / n), k) + fx(a + (b - a) * (i / n), k)) / (2 * n)
        return 1 - summa

    d, nj = intervals(K)
    qj = q(d)
    R0 = R0(qj, nj)
    F_ = F_(R0, K)
    if (F_ > alp):
        print('Принимается')
    else:
        print('Не принимается')


def third():
    ncount1 = []
    ncount2 = []

    for i in range(100):
        r = int(r_value.get())
        N1 = int(N_value.get())
        ver1 = float(ver1_value.get())
        ver2 = float(ver2_value.get())
        ver3 = float(ver3_value.get())
        K = int(k_value.get())
        alp = float(alpha_value.get())
        arr_sum = []

        for i in range(0, N1):
            sum1 = 1
            sum2 = 1
            sum3 = 1
            flag1 = 0
            flag2 = 0
            flag3 = 0
            for j in range(0, r - 1):
                p1 = random.random()
                p2 = random.random()
                p3 = random.random()
                if (p1 > ver1 and flag1 == 0):
                    sum1 = sum1 + 1
                else:
                    flag1 = 1
                if (p2 > ver2 and flag2 == 0):
                    sum2 = sum2 + 1
                else:
                    flag2 = 1
                if (p3 > ver3 and flag3 == 0):
                    sum3 = sum3 + 1
                else:
                    flag3 = 1
            sum = sum1 + sum2 + sum3
            arr_sum.insert(i, sum)

        c = Counter(arr_sum)
        c.most_common()
        c = sorted(c.items())
        arr_sum.sort()
        Arr_sum = []
        for i in range(r, 3*r+1):
            Arr_sum.append(i)
        lst = []
        for i in range(len(c)):
            lst.append(c[i][0])
        ni = []
        j = 0
        if lst[0] != Arr_sum[0]:
            ni.append(0)
            j += 1
        i = 0
        while len(lst) != len(Arr_sum):
            lst.append(0)
        while j < len(Arr_sum) and i < len(lst):
            if (Arr_sum[j] == lst[i]):
                ni.append(c[i][1])
                i += 1
                j += 1
            else:
                ni.append(0)
                j += 1
        ver1_ = []
        ver2_ = []
        ver3_ = []

        for i in range(1, r + 1):
            ver_1 = ((1 - ver1) ** (i - 1)) * ver1
            ver_2 = ((1 - ver2) ** (i - 1)) * ver2
            ver_3 = ((1 - ver3) ** (i - 1)) * ver3
            if (i == r):
                ver_1 += (1 - ver1) ** r
                ver_2 += (1 - ver2) ** r
                ver_3 += (1 - ver3) ** r
            ver1_.append(ver_1)
            ver2_.append(ver_2)
            ver3_.append(ver_3)

        vers = []
        for k in range(3, 3 * r + 1):
            vers_ = 0
            for i in range(1, r + 1):
                for j in range(1, r + 1):
                    if ((k - i - j - 1) >= 0 and (k - i - j) <= len(ver3_)):
                        vers_ += ver1_[i - 1] * ver2_[j - 1] * ver3_[k - i - j - 1]
            vers.append(round(vers_, 9))
        def Gamma_func(value):
            if (value == 1):
                return value
            if (value == 0.5):
                return np.sqrt(np.pi)
            return ((value - 1) * Gamma_func(value - 1))


        def intervals(k, ounter_add):
            tmp = vers
            k+=1
            num = len(tmp) / k
            t = []
            d = []
            ntmp = []
            nj_ = []
            nj = []
            count = 0
            count_2 = num
            sum = 0
            for i in tmp:
                if (ounter_add == k - 1):
                    t.append(i)
                    ntmp.append(ni[count])
                    count += 1
                    sum += i
                    continue
                t.append(i)
                ntmp.append(ni[count])
                count += 1
                sum+=i
                #print(i, 'i')
                print(sum)

                if (sum >= 1/k):
                    d.append(t)
                    nj_.append(ntmp)
                    ntmp = []
                    t = []
                    sum=0
                    count_2 += num
                    ounter_add +=1

            d.append(t)
            nj_.append(ntmp)

            for i in nj_:
                nj.append(np.sum(i))
            return d, nj

        def intervals1(k):
            tmp = vers
            num = len(tmp) / k
            t = []
            d = []
            ntmp = []
            nj_ = []
            nj = []
            count = 0
            count_2 = num
            # print(num)
            for i in tmp:
                t.append(i)
                ntmp.append(ni[count])
                count += 1
                if (count >= count_2):
                    # print(count)
                    # print(count_2)
                    # print(t)
                    d.append(t)
                    nj_.append(ntmp)
                    ntmp = []
                    t = []
                    count_2 += num
            # print("---",nj_)
            for i in nj_:
                nj.append(np.sum(i))
            return d, nj

        def q(d):
            q = []
            for i in d:
                q.append(np.sum(i))
            return q

        def R0(qj, nj):
            summa = 0
            for i in range(len(qj)):
                summa += (nj[i] - N1 * qj[i]) ** 2 / (N1 * qj[i])
            return summa

        def fx(x, r):
            if (x > 0):
                return (2 ** (-r / 2)) * (Gamma_func(r / 2) ** (-1)) * (x ** (r / 2 - 1) * (np.exp(-x / 2)))
            else:
                return 0

        def F_(R0, k, n=1000):
            a = 0
            b = R0
            summa = 0
            k = k - 1
            for i in range(1, n + 1, 1):
                summa += (b - a) * (fx(a + (b - a) * ((i - 1) / n), k) + fx(a + (b - a) * (i / n), k)) / (2 * n)
            return 1 - summa

        d, nj = intervals1(K)
        qj = q(d)
        R0 = R0(qj, nj)
        F_ = F_(R0, K)
        if (F_ > alp):
            ncount1.append(1)
        else:
            ncount2.append(1)
    print(qj)
    pr = Label(text='Принялось: ' + str(np.sum(ncount1)))
    pr.grid(row=13, column=0, columnspan=2, padx=10, pady=0, sticky="we")
    ot = Label(text='Отверглось: ' + str(np.sum(ncount2)))
    ot.grid(row=14, column=0, columnspan=2, padx=10, pady=0, sticky="we")

def make_table():
    Names = ['y', 'ni', 'ni/n']
    N_ = int(N_value.get())

    table.grid(row=0, column=3, rowspan=7, columnspan=1, padx=10, pady=10, sticky="ew")
    table["columns"] = Names
    for head in Names:
        table.heading(head, text=head, anchor=CENTER)
        table.column(head, anchor=CENTER)
        table.column(head, width=250)

    scroll_bar1 = Scrollbar(window, orient='vertical', command=table.yview)
    scroll_bar1.grid(row=0, column=4, rowspan=7, columnspan=1, pady=10, sticky=NSEW)
    table.configure(yscroll=scroll_bar1.set)

def make_table_2():
    Names1 = ['E' + chr(951), 'x\u0304', '|' + 'E' + chr(951) + ' - ' + 'x\u0304' + '|', 'D' + chr(951), 'S\u00B2', '|' + 'D' + chr(951) + ' - ' + 'S\u00B2' + '|', 'M\u0302' + 'e', 'R\u0302']
    N_ = int(N_value.get())

    table1.grid(row=7, column=3, rowspan=7, columnspan=1, padx=10, pady=10, sticky="nsew")
    table1["columns"] = Names1
    for head in Names1:
        table1.heading(head, text=head, anchor=CENTER)
        table1.column(head, anchor=CENTER)
        table1.column(head, width=25)

    '''scroll_bar2 = Scrollbar(window, orient='vertical', command=table1.yview)
    scroll_bar2.grid(row=7, column=4, rowspan=7, columnspan=1, pady=10, sticky=NSEW)
    table1.configure(yscroll=scroll_bar2.set)'''

def make_table_3():
    Names2 = ['y', 'P(y)', 'nj / n']
    N_ = int(N_value.get())

    table2.grid(row=15, column=3, rowspan=7, columnspan=1, padx=10, pady=10, sticky="nsew")
    table2["columns"] = Names2
    for head in Names2:
        table2.heading(head, text=head, anchor=CENTER)
        table2.column(head, anchor=CENTER)
        table2.column(head, width=25)

    scroll_bar2 = Scrollbar(window, orient='vertical', command=table2.yview)
    scroll_bar2.grid(row=15, column=4, rowspan=7, columnspan=1, pady=10, sticky=NSEW)
    table.configure(yscroll=scroll_bar2.set)

def refresh_table1(c):
    N_ = int(N_value.get())
    for row in table.get_children():
        table.delete(row)
    for i in range(len(list(c))):
        table.insert('', END, values=(int(c[i][0]), int(c[i][1]), float(c[i][1]/N_)))


def refresh_table2(c):
    N_ = int(N_value.get())
    for row in table1.get_children():
        table1.delete(row)
    table1.insert('', END, values=(round(c[0], 7), round(c[1], 7), round(c[2], 7), round(c[3], 7), round(c[4], 7), round(c[5], 7), round(c[6], 7), round(c[7], 7)))

def refresh_table3(c, vers):

    N_ = int(N_value.get())
    for row in table2.get_children():
        table2.delete(row)
    for i in range(len(list(c))):
        table2.insert('', END, values=(int(c[i][0]), vers[i], float(c[i][1] / N_)))


r_name = StringVar()
N_name = StringVar()
ver1_name = StringVar()
ver2_name = StringVar()
ver3_name = StringVar()
N2_name = StringVar()
k_name = StringVar()
alpha_name = StringVar()


N_ = Label(text="Введите количество экспериментов: ")
N_.grid(row=0, column=0, sticky="w", padx=10, pady=10)
N_value = Entry(textvariable=N_name)
N_value.grid(row=0, column=1, padx=5, pady=5)
N_value.insert(0, "500")

r = Label(text="Введите количество выстрелов: ")
r.grid(row=1, column=0, sticky="w", padx=10, pady=10)
r_value = Entry(textvariable=r_name)
r_value.grid(row=1, column=1, padx=5, pady=5)
r_value.insert(0, "3")

ver1 = Label(text="Введите вер-ть попадания 1 стрелка: ")
ver1.grid(row=2, column=0, sticky="w", padx=10, pady=10)
ver1_value = Entry(textvariable=ver1_name)
ver1_value.grid(row=2, column=1, padx=5, pady=5)
ver1_value.insert(0, "0.8")

ver2 = Label(text="Введите вер-ть попадания 2 стрелка: ")
ver2.grid(row=3, column=0, sticky="w", padx=10, pady=10)
ver2_value = Entry(textvariable=ver2_name)
ver2_value.grid(row=3, column=1, padx=5, pady=5)
ver2_value.insert(0, "0.8")

ver3 = Label(text="Введите вер-ть попадания 3 стрелка: ")
ver3.grid(row=4, column=0, sticky="w", padx=10, pady=10)
ver3_value = Entry(textvariable=ver3_name)
ver3_value.grid(row=4, column=1, padx=5, pady=5)
ver3_value.insert(0, "0.8")

k = Label(text="Введите число интервалов: ")
k.grid(row=5, column=0, sticky="w", padx=10, pady=10)
k_value = Entry(textvariable=k_name)
k_value.grid(row=5, column=1, padx=5, pady=5)
k_value.insert(0, "3")

alpha = Label(text="Введите уровень значимости: ")
alpha.grid(row=6, column=0, sticky="w", padx=10, pady=10)
alpha_value = Entry(textvariable=alpha_name)
alpha_value.grid(row=6, column=1, padx=5, pady=5)
alpha_value.insert(0, "0.5")

message_button = Button(text="Решить", command=run)
message_button.grid(row=7, column=0, columnspan=2, rowspan=1, padx=10, pady=5, sticky="nsew")


task = Button(text='Задача', bg='#ececec', highlightbackground='#ececec', command=task_window, width = 15)
task.grid(row=8, column=0, columnspan=2, padx=10, pady=5, sticky="we")

graphs = Button(text='Графики', bg='#ececec', highlightbackground='#ececec', command=graphs_window, width = 15)
graphs.grid(row=9, column=0, columnspan=2, padx=10, pady=0, sticky="we")

tr = Button(text='Гипотеза', bg='#ececec', highlightbackground='#ececec', command=third, width = 15)
tr.grid(row=12, column=0, columnspan=2, padx=10, pady=0, sticky="we")

table = ttk.Treeview(window, show='headings', selectmode='browse')
table1 = ttk.Treeview(window, show='headings', selectmode='browse')
table2 = ttk.Treeview(window, show='headings', selectmode='browse')

make_table()
make_table_2()
make_table_3()
window.mainloop()
