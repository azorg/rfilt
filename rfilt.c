/*
 * Проект: "Нелинейный фильтр ограничения скорости, ускорения и рывка"
 * Версия: 0.1a
 * Файл: "rfilt.c"
 */

//-----------------------------------------------------------------------------
#include <math.h>
#include "rfilt.h"
//-----------------------------------------------------------------------------
// инициализация фильтра
void rfilt_init(
  rfilt_t *self,     // указатель на внутреннюю структуру фильтра
  double dt,         // шаг вызова процедуры фильтрации по времени
  double v_max,      // максимальная скорость
  double a_up_max,   // максимальное ускорение разгона
  double a_down_max, // максимальное ускорение торможения
  double r_max,      // максимальный рывок
  double x_init,     // начальное положение
  double x_min,      // ораничение снизу
  double x_max)      // ораничение сверху
{
  rfilt_tune(self, dt, v_max, a_up_max, a_down_max, r_max, x_min, x_max);

  self->x = x_init;
  self->v = 0.;
  self->a = 0.;

  self->s  = 0.;
  self->nt = 0.;
}
//-----------------------------------------------------------------------------
// установка параметров фильтра "на лету"
void rfilt_tune(
  rfilt_t *self,     // указатель на внутреннюю структуру фильтра
  double dt,         // шаг вызова процедуры фильтрации по времени
  double v_max,      // максимальная скорость
  double a_up_max,   // максимальное ускорение разгоа
  double a_down_max, // максимальное ускорение разгоа
  double r_max,      // максимальный рывок
  double x_min,      // ораничение снизу
  double x_max)      // ораничение сверху
{
  double dt2 = dt * dt;
  self->V    = v_max * dt;
  self->Au   = a_up_max   * dt2;
  self->Ad   = a_down_max * dt2;
  self->R    = r_max * dt2 * dt;
  self->Xmin = x_min;
  self->Xmax = x_max;
}
//-----------------------------------------------------------------------------
// выполнить шаг фильтрации
// (функция возвращает "профильтрованное" значение с учетом ограничений)
double rfilt_step(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра
  double x)      // входное значение для фильтра
{
  double trend, A, ds, r, v, a;
  double rmin = -self->R, rmax = self->R;
  double a2 = self->a * self->a;
  double vcrit, s, nt, smin, smax;
  
  // ограничить заданное значение
  x = RFILT_LIMIT(x, self->Xmin, self->Xmax);
  
  trend = self->v * self->a;
  if (trend < 0.)
  { // торможение (скорость и ускорение не равны нулю и имеют разный знак
    A = self->Ad; // максимальное приведенное ускорение торможения
  }
  else // if (trend >= 0)
  { // разгон (знаки скорости и ускорения совпадают)
    A = self->Au; // максмальное приведенное ускорение разгона
  }

  // вычислить тормозной путь и время тормозного пути
  s = rfilt_s(self, self->v, self->a, &nt);
  self->s  = s;
  self->nt = nt;

  // требуемое перемещение
  ds = x - self->x;

#if 0 // FIXME
  // проверить возможность "малого хода"
  r = (ds - self->v - self->a * 0.5) * 6.; // рывок для перемещения за 1 такт
  v = self->v + self->a + r * 0.5; // скорость в конце такта
  a = self->a + r; // ускорение в конце такта
  
  if (fabs(r) <= self->R && fabs(v) <= self->V && fabs(a) <= A)  
  { // есть условия для малого хода!
#if 0 // педантично
    self->x += self->v + self->a * 0.5 + r / 6.;
    self->v += self->a + r * 0.5;
    self->a += r;
#else // упрощенно
    self->x = x;
    self->v = v;
    self->a = a;
#endif
    RFILT_DBG("pipe mode: x=%f v=%f a=%f r=%f", self->x, self->v, self->a, r);
    return self->x;
  }
#endif

  // установить ограничения на величину рывка
  // огриничить рывок в связи с ограничением ускорения
  if (rmax >  A - self->a) rmax =  A - self->a;
  if (rmin < -A - self->a) rmin = -A - self->a;
  
#if 1 // FIXME 
  // вычислить критическую скорость, выше которой требуется
  // "быстрое" снижения ускорения разгона для предотращения
  // превышения по модулю максимальной скорости self->V
  vcrit = self->V - a2 / (2. * self->R);
  
  if ((self->a >= 0. && self->v >=  vcrit) || // скорость растёт и достигнет V
      (self->a <= 0. && self->v <= -vcrit))   // скорость падает и достигнет -V
  { // достигнута критическая скорость, требуется немедленное обнуление ускорения
    // во избежание превышения максимальной скорости
    // задаем максимальный рывок с целью снижения ускорения по модулю
    if (self->a != 0.)
    {
      r = RFILT_LIMIT(-self->a, rmin, rmax);
      self->x += self->v + self->a * 0.5 + r / 6.;
      self->v += self->a + r * 0.5;
      self->a += r;
      RFILT_DBG("limit Vmax: x=%f v=%f a=%f r=%f", self->x, self->v, self->a, r);
      return self->x;
    }
    else
    { // если критеческая скорость достигнута, а=0 установить запрет на рывок
      if (self->v > 0.) rmax = 0.;
      if (self->v < 0.) rmin = 0.;
    }
  }
#endif

#if 1 // FIXME
  // проверка возможности торможения до ограничителей Xmax и Xmin
  if (self->x + s >= self->Xmax)
  {
    r = rmin;
    self->x += self->v + self->a * 0.5 + r / 6.;
    self->v += self->a + r * 0.5;
    self->a += r;
    RFILT_DBG("stop over Xmax (!): x=%f v=%f a=%f r=%f", self->x, self->v, self->a, r);
    return self->x;
  }
  else if (self->x + s <= self->Xmin)
  {
    r = rmax;
    self->x += self->v + self->a * 0.5 + r / 6.;
    self->v += self->a + r * 0.5;
    self->a += r;
    RFILT_DBG("stop over Xmin (!): x=%f v=%f a=%f r=%f", self->x, self->v, self->a, r);
    return self->x;
  }
#endif

#if 0 // FIXME
  // проверить возможность "разгона" и пропорционального управления рывком
  // если (x[i]+s) "перепрыгнул" через заданный x => тормозим без аналитики
  if ((ds >= 0. && s >= ds) ||
      (ds  < 0. && s <= ds))
  {
    if (s >= 0.) r = rmin;
    else         r = rmax;
    self->x += self->v + self->a * 0.5 + r / 6.;
    self->v += self->a + r * 0.5;
    self->a += r;
    RFILT_DBG("stop over: x=%f v=%f a=%f r=%f", self->x, self->v, self->a, r);
    return self->x;
  }
#endif

  // пропорциональное управление рывком
  v = self->v + self->a + rmin * 0.5;
  a = self->a + rmin;
  smin = self->v + self->a * 0.5 + rmin / 6. + rfilt_s(self, v, a, &nt);
  v = self->v + self->a + rmax * 0.5;
  a = self->a + rmax;
  smax = self->v + self->a * 0.5 + rmax / 6. + rfilt_s(self, v, a, &nt);
  if ((ds > smin && ds < smax) ||
      (ds > smax && ds < smin))
  { // пропорциональное управление рывком
    r = (rmax - rmin) / (smax - smin) * (ds - smin) + rmin;
  }
  else
  { // не пропорциональный разгон/торможение
    if      (ds >= smin && ds >= smax) r = rmax;
    else if (ds <= smin && ds <= smax) r = rmin;
    else                               r = 0.;
  }

  self->x += self->v + self->a * 0.5 + r / 6.;
  self->v += self->a + r * 0.5;
  self->a += r;

  RFILT_DBG("work: x=%f v=%f a=%f r=%f s=%f smin=%f smax=%f rmin=%f rmax=%f",
             self->x, self->v, self->a, r, s, smin, smax, rmin, rmax);
  return self->x;
}
//-----------------------------------------------------------------------------
// расчёт тормозного пути
// (функция возвращяет минимальный тормозной путь со знаком, а так же nt)
double rfilt_s(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double *nt)    // минимальное время торможения в тактах
{
  double s = 0., a2 = a * a, Ad2 = self->Ad * self->Ad;
  double dt, dt2;
  //double R2 = self->R * self->R, Ad3 = Ad2 * self->Ad;
  double sig = (v >= 0.) ? 1. : -1.;
  v *= sig; // теперь v >= 0
  a *= sig;

  if (a >= 0.)
  { // тормозной путь после разгона (a и v одного знака) или
    // после равномерного движения с постоянной скоростью v и ускорением a=0.
    
    // вычислить максимальную приведенную скорость (фактически по модулю)
    // достигаемую за время обнуления начального укорения разгона
    // с максимальным рывком
    double vmax = v + a2 / (2. * self->R); // vmax >= 0

    dt  = a / self->R; // время обнуления текущего ускорения
    s   = dt * (v + a2 / (3. * self->R)); // тормозной путь обнуления ускорения
    *nt = dt;

    if (vmax <= self->Ad * self->Ad / self->R)
    { // торможение состоит из 2-х фаз
      dt   = sqrt(vmax / self->R); // время каждой фазы
      s   += dt * vmax; // тормозной путь ф1 и ф2 
      *nt += dt * 2.;   // время торможения ф1 и ф2
    }
    else
    { // торможение состоит из 3-х фаз
      dt   = vmax / self->Ad + self->Ad / self->R; // время всех фаз
      s   += vmax * dt / 2.; // тормозной путь ф1..ф3
      *nt += dt;             // время торможения ф1..ф3
    }
  }
  else // if (a < 0.)
  { // тормозной путь после начала торможения
    
    // вычислить критические пороги по скорости
    double v1 = a2 / (self->R * 2.); // с меньшей скоростю торможение с выбегом 
    double v2 = Ad2 / self->R - v1;

    if (v >= v2)
    { // торможение в 3 фазы
      
      // фаза 1 (увеличиваем тормозное ускорения с а до Ad)
      dt  = (self->Ad + a) / self->R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. - self->R * dt2 * dt / 6.;
      *nt = dt; 

      // фаза 2 (тормозим с постоянным ускорением Ad)
      dt   = (v - v2) / self->Ad;
      s   += (v + a2 / (2. * self->R)) * dt * 0.5;
      *nt += dt;

      // фаза 3 (обнуляем скорость и ускорение с максимальным рывком)
      dt   = self->Ad / self->R;
      s   += self->R * dt * dt * dt / 6.;
      *nt += dt;
    }
    else if (v >= v1) // && v < v2
    { // торможение в 2 фазы

      // максимальное по модулу ускорение в конце первой фазы:
      double amax = sqrt(self->R * v + a2 * 0.5);
      
      // фаза 1 (увеличиваем тормозное ускорения с а до amax)
      dt = (amax + a) / self->R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. - self->R * dt2 * dt / 6.;
      *nt = dt; 
      
      // фаза 2 (обнуляем скорость и ускорение с максимальным рывком)
      dt   = amax / self->R;
      s   += self->R * dt * dt * dt / 6.;
      *nt += dt;
    }
    else // if (v < v1)
    { // торможение с "выбегом" (3 или 4 фазы)

      // фаза 1 (сбрасываем ускорение до нуля)
      dt = -a / self->R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. + self->R * dt2 * dt / 6.;
      *nt = dt; 
      
      // теперь задача сводится к торможению после равномерного движения
      // в две или три фазы
      v = v1 - v; // считаем, что остаточная скорость положительна для удобства

      if (v <= Ad2 / self->R)
      { // торможение состоит из 3-х фаз
        dt   = sqrt(v / self->R); // время фазы 2 или 3
        s   -= dt * v;  // тормозной путь ф2 и ф3 
        *nt += dt * 2.; // время торможения ф2 и ф3
      }
      else
      { // торможение состоит из 4-х фаз
        dt   = v / self->Ad + self->Ad / self->R; // время фаз 2, 3 и 4
        s   -= v * dt / 2.; // тормозной путь ф2..ф4
        *nt += dt;          // время торможения ф2..ф4
      }
    }
  }

  return s * sig; // вернуть тормозной путь со знаком
}
//-----------------------------------------------------------------------------

/*** end of "rfilt.c" file ***/


