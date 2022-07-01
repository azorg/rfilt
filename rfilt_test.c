/*
 * Проект: "Нелинейный фильтр ограничения скорости, ускорения и рывка"
 * Файл: "rfilt_test.h" (тестовый проект)
 */

//-----------------------------------------------------------------------------
#include <stdlib.h> // malloc(), free()
#include <stdio.h>  // fopen(), fprintf(), fclose()
#include <math.h>   // M_PI, sin()
#include "rfilt.h"
//-----------------------------------------------------------------------------
int main()
{
  int i, n = 4000;   // число точек
  double dt = 0.010; // 10 мс - шаг разбиения
  FILE *f = fopen("data.txt", "wb"); // файл анных

  double xp, vp, ap, *d = malloc(sizeof(double) * n);
  rfilt_t flt;

  // инициализировать фильтр
  rfilt_init(
    &flt,    // указатель на внутреннюю структуру фильтра
    dt,      // шаг вызова процедуры фильтрации по времени
    70.,     // максимальная скорость
    80.,     // максимальное ускорение разгона
    85.,     // максимальное ускорение торможения
    120.,    // максимальный рывок
    0.,      // коэфицент бета фильтра [0...1]
    0.,      // начальное положение
    -600.,   // ораничение снизу
    600.);   // ораничение сверху

  // заполнить массив входных значений
  for (i = 0; i < n; i++) d[i] = 0.;
  for (i = 20; i < 820; i++) d[i] = 150.;

#if 1
  for (i = 820; i < 1120; i++) d[i] = -30.;

  for (i = 1400; i < (1400+1000); i++) // 2 периода
    d[i] = 20. * sin((i - 1400.) * 2 * M_PI / 500.); // период 5 сек

  for (i = 2800; i < 2900; i++) d[i] = i - 2800;
  for (i = 2900; i < 3000; i++) d[i] = 3000 - i;

  for (i = 3500; i < 4000; i++) d[i] = -100;
#endif

  // прогнать фильтр
  xp = vp = ap = 0.;
  for (i = 0; i < n; i++)
  //for (i = 0; i < 1000; i++)
  {
    double y = d[i];
    double x = rfilt_step(&flt, y);
    double v = (x - xp) / dt;
    double a = (v - vp) / dt;
    double r = (a - ap) / dt;
    //printf("<<< %8i\n", i);

    fprintf(f, "%8i %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n",
            i, y, x, v, a, r, RFILT_LIMIT_ABS(xp + flt.s, 10000.));
    xp = x;
    vp = v;
    ap = a;
  }

  free(d);

  return 0;
}
//-----------------------------------------------------------------------------

/*** end of "rfilt_test.c" file ***/

