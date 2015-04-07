# Алгоритм

0. Используем исходное пространственное положение последовательности аминокислот, рассматриваем пару цепочек.
Считаем, что входные данные приведены в PDB-файле.
1. Определяем модель, описывающую относительные положения элементов цепочек в пространстве.
2. Выбираем критерий, по которому будем отбирать аминокислоты для мутагенеза.
- расстояние Хаусдорфа
- ближайший атом
- что-то с узлами (что?)
3. Отобранные на основании критерия аминокислоты подвергаем мутагенезу: заменяем на аланин, если электростатические свойства взаимодействующих цепочек улучшились, то подбираем вместо аланина ту аминокислоту, которая улучшает эти свойства сильнее всего.

# Детализация (что есть сейчас в том или ином виде)

1. Модель

- аминокислота представлена координатами всех атомов
- строим по этим данным трехмерную триангуляцию Делоне, в которой по каждому атому можно восстановить аминокислоту

- для построения триангуляции Делоне используется библиотека qhull.

! подумать, может быть стоит использовать какое-то упрощение модели? Если нет, то почему

2. Критерий

  2.1. Расстояние Хаусдорфа, или поиск ближайшей пары атомов
  Идея использовать расстояние Хаусдорфа в качестве метрики для определения ближайших аминокислот - из статьи   (вспомнить, откуда):
  цепочки представлены как поверхности, проходящие через центры атомов. Расстояние определяется между поверхностями и в терминах поверхностей - отсюда идея использовать метрику, определяющую пространственное расстояние.
  Подход позволяет учесть пространственные свойства атомов: случай взаимодействия атома и треугольной грани, случай взаимодействия двух пар атомов.  
  2.2.

3. Мутагенез: монте-карло

- для мутагенеза используется классическая версия алгоритма монте-карло, сейчас используется готовая реализация протокола  из pyrosetta.