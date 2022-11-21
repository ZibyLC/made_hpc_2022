**Знакомство в компании**

Представим, что процессы – это компания незнакомых людей, которые знакомятся с помощью следующей игры:
<ol>
<li> Начинает процессор 0. Случайным образом он выбирает другой процессор i и посылает ему сообщение со своим именем (можете случайным образом задавать имя)
<li> Процессор i отсылает сообщение случайному процессору j (которые еще не участвовал в игре), в сообщении – все имена и ранги предыдущих процессоров в правильном порядке. Номер процессора j знает только I, так что все должны быть начеку.
<li> Игра заканчивается через N ходов. Используйте синхронную пересылку MPI_SSend
</ol>
Напишите программу используя MPI. (25 баллов) - **Решение в папке acquaintance**


**Параллельный одномерный клеточный автомат.**

С помощью MPI распараллельте одномерный клеточный автомат Вольфрама (Rule110).
Игра происходит следующим образом:
<ol>
<li> Инициализируйте одномерный массив 0 и 1 случайным образом
<li> В зависимости от значений: левого соседа, себя, правого соседа на следующем шаге клетка либо меняет значение, либо остается той же. Посмотрите, например, что значит Rule110 (https://en.wikipedia.org/wiki/Rule_110)
</ol>
Сделайте периодические и непериодические граничные условия (5 баллов)
Работает параллельный код на нескольких процессах (20 баллов)
Имплементированы клетки-призраки (ghost cells) (10 балла)
Можно поменять правило игры (сделать одно из 256) (20 баллов)
График ускорения работы программы от кол-ва процессов (5 баллов)
Картинка эволюции для одного правила (15 баллов)
**Решение в папке automaton**
Итого баллов: 75  + 25 = 100 баллов за базовую часть из 2 заданий