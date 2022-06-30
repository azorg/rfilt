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
  double beta,       // коэфицент бета фильтра [0...1]
  double x_init,     // начальное положение
  double x_min,      // ораничение снизу
  double x_max)      // ораничение сверху
{
  rfilt_tune(self, dt, v_max, a_up_max, a_down_max, r_max, beta, x_min, x_max);

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
  double beta,       // коэфицент бета фильтра [0...1]
  double x_min,      // ораничение снизу
  double x_max)      // ораничение сверху
{
  double dt2 = dt * dt;
  self->V    = v_max * dt;
  self->Au   = a_up_max   * dt2;
  self->Ad   = a_down_max * dt2;
  self->R    = r_max * dt2 * dt;
  self->B    = RFILT_LIMIT(beta, 0., 1.);
  self->Xmin = x_min;
  self->Xmax = x_max;
}
//-----------------------------------------------------------------------------
// ограничение рывка в связи с ограничением скорости и ускорния
static double rfilt_limit_r(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double r)      // предлагаемый рывок
{
  double A, vc;

  if (v * a < 0.)
  { // торможение (скорость и ускорение не равны нулю и имеют разный знак
    A = self->Ad; // максимальное приведенное ускорение торможения
  }
  else // if (trend >= 0)
  { // разгон (знаки скорости и ускорения совпадают или нули)
    A = self->Au; // максмальное приведенное ускорение разгона
  }

#if 1 // ограничить рывок в связи с ограничением ускорения
  if (a + r >  A) r =  A - a;
  if (a + r < -A) r = -A - a;
  // возмозна ситуация смены знака скорости при этом возникает
  // большой рывок с Ad на Au (или наоборот)
  //r = RFILT_LIMIT_ABS(r, self->R); // FIXME: !!!???
#endif

#if 1 // ограничить рывок в связи с ограничением скорости
  // вычислить критическую скорость, выше которой требуется
  // "быстрое" снижения ускорения разгона для предотращения
  // превышения по модулю максимальной скорости self->V
  vc  = self->V - a  * a  / (2. * self->R);

  if ((a > 0. && v >=  vc) ||
      (a < 0. && v <= -vc))
  { // текущая скорость достигла критической
    r = RFILT_LIMIT_ABS(-a, self->R);
    RFILT_DBG("critical V: v=%f a=%f r=%f vcrit=%f vapprox=%f",
              v, a, r, vc, v + a * a / (2. * self->R));
  }
  else if (a >= 0. && v > self->V)
  { // скорость превысила +Vmax
    a = 0.;
    if (r > 0.) r = 0.;
    RFILT_DBG("v=%f > Vmax=%f r=%f", v, self->V, r);
  }
  else if (a <= 0. && v < -self->V)
  { // скорость ниже -Vmin
    a = 0.;
    if (r < 0.) r = 0.;
    RFILT_DBG("v=%f < -Vmax=%f r=%f", v, -self->V, r);
  }
#endif

  return r;
}
//----------------------------------------------------------------------------
// расчёт тормозного пути
// (функция возвращяет минимальный тормозной путь со знаком, а так же nt)
static double rfilt_s(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double R,      // приведенное ограничение по рывку
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double *nt)    // минимальное время торможения в тактах
{
  double s = 0., a2 = a * a, Ad2 = self->Ad * self->Ad;
  double dt, dt2;
  double sig = (v >= 0.) ? 1. : -1.;
  v *= sig; // теперь v >= 0
  a *= sig;

  if (a >= 0.)
  { // тормозной путь после разгона (a и v одного знака) или
    // после равномерного движения с постоянной скоростью v и ускорением a=0.

    // вычислить максимальную приведенную скорость (фактически по модулю)
    // достигаемую за время обнуления начального укорения разгона
    // с максимальным рывком
    double vmax = v + a2 / (2. * R); // vmax >= 0

    dt  = a / R; // время обнуления текущего ускорения
    s   = dt * (v + a2 / (3. * R)); // тормозной путь обнуления ускорения
    *nt = dt;

    if (vmax <= self->Ad * self->Ad / R)
    { // торможение состоит из 2-х фаз
      dt   = sqrt(vmax / R); // время каждой фазы
      s   += dt * vmax; // тормозной путь ф1 и ф2
      *nt += dt * 2.;   // время торможения ф1 и ф2
    }
    else
    { // торможение состоит из 3-х фаз
      dt   = vmax / self->Ad + self->Ad / R; // время всех фаз
      s   += vmax * dt / 2.; // тормозной путь ф1..ф3
      *nt += dt;             // время торможения ф1..ф3
    }
  }
  else // if (a < 0.)
  { // тормозной путь после начала торможения

    // вычислить критические пороги по скорости
    double v1 = a2 / (R * 2.); // с меньшей скоростю торможение с выбегом
    double v2 = Ad2 / R - v1;

    if (v >= v2)
    { // торможение в 3 фазы c набором тормозного ускорения Ad

      // фаза 1 (увеличиваем тормозное ускорения с а до Ad)
      dt  = (self->Ad + a) / R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. - R * dt2 * dt / 6.;
      *nt = dt;

      // фаза 2 (тормозим с постоянным ускорением Ad)
      dt   = (v - v2) / self->Ad;
      s   += (v + a2 / (2. * R)) * dt * 0.5;
      *nt += dt;

      // фаза 3 (обнуляем скорость и ускорение с максимальным рывком)
      dt   = self->Ad / R;
      s   += R * dt * dt * dt / 6.;
      *nt += dt;
    }
    else if (v >= v1) // && v < v2
    { // торможение в 2 фазы c |amax|<Ad

      // максимальное по модулу ускорение в конце первой фазы
      double amax = sqrt(R * v + a2 * 0.5);

      // фаза 1 (увеличиваем тормозное ускорения с а до amax)
      dt = (amax + a) / R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. - R * dt2 * dt / 6.;
      *nt = dt;

      // фаза 2 (обнуляем скорость и ускорение с максимальным рывком)
      dt   = amax / R;
      s   += R * dt * dt * dt / 6.;
      *nt += dt;
    }
    else // if (v < v1)
    { // торможение с "выбегом" (3 или 4 фазы)

      // фаза 1 (сбрасываем ускорение до нуля)
      dt = -a / R; // помним, что a<0
      dt2 = dt * dt;
      s   = v * dt + a * dt2 / 2. + R * dt2 * dt / 6.;
      *nt = dt;

      // теперь задача сводится к торможению после равномерного движения
      // в две или три фазы
      v = v1 - v; // считаем, что остаточная скорость положительна для удобства

      if (v <= Ad2 / R)
      { // торможение состоит из 3-х фаз
        dt   = sqrt(v / R); // время фазы 2 или 3
        s   -= dt * v;  // тормозной путь ф2 и ф3
        *nt += dt * 2.; // время торможения ф2 и ф3
      }
      else
      { // торможение состоит из 4-х фаз
        dt   = v / self->Ad + self->Ad / R; // время фаз 2, 3 и 4
        s   -= v * dt / 2.; // тормозной путь ф2..ф4
        *nt += dt;          // время торможения ф2..ф4
      }
    }
  }

  return s * sig; // вернуть тормозной путь со знаком
}
//-----------------------------------------------------------------------------
#if 0
// расчёт рывка
double rfilt_r_bad(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double dx)     // невязка заданой и текущей позиции
{
  double r, s;
  double R = self->R;
  double sig = (dx >= 0.) ? 1. : -1.;

  // вычислить тормозной путь и время тормозного пути с максимальным рывком
  s = self->s = rfilt_s(self, R, v, a, &self->nt);

  // привсти знаки относитльн dx
  v  *= sig;
  a  *= sig;
  s  *= sig;
  dx *= sig; // теперь dx >= 0

  // пропорциональное управление рывком
  //R *= RFILT_LIMIT(fabs(dx) / 10., 0., 1.); // !!! FIXME: magic

  if (s < 0.)
  { // s < x <= dx - "реверс"
    // v>0 a<0 - неверное торможение
    // v<0 a>0 - вход в реверс
    // v<0 a<0 - провал в реверс
    // v>0 a>0 - невероятный разгон (не может быть!)
    r = R;
    RFILT_DBG("reverse s<0");
  }
  else if (s > dx)
  { // s <= dx <= s -"пролёт"
    // v>0 a<0 - торможение
    // v<0 a>0 - выход из реверса
    // v<0 a<0 - провал в реврс
    // v>0 a>0 - разгон
#if 1 // FIXME: !!!
    // торможение состоит из от 1-й до 4-х фаз:
    // 1. сброс ускорения разгона с рывком r<0 и набором сокрости (если a>0)
    // 2. набор (добор) тормозного ускорения со сбросом скорости и с рывком r<0
    // 3. сброс скорости с постоянным ускорением без рывка r=0 (если ускорение достигло максимума)
    // 4. сброс скорости и тормозного ускорения до нуля с положительным рывком r>0
    double a2 = a * a;
    double Ad2 = self->Ad * self->Ad;
    if (a > 0.)
    { // фаза 1: требуется сброс ускорения разгона до 0 и далее
      r = -R;
      RFILT_DBG("stop phase 1");
    }
    else // if a <= 0
    {
      if (2. * v * R <= a2)
      { // фаза 4: требуется сброс тормозного ускорения и остаточной скорости до нуля
        //r = R;
        r = a2 / (2. * v); // FIXE
        if (r > -a) r = a;
        // FIXME:
        //!!!self->v = a2 / (2. * R); // коррекция накапливаемых ошибок вычислений!
        RFILT_DBG("stop phase 4");
      }
      else if (v <= v + (a2 - Ad2) / (2. * R) ||
               a <= -self->Ad)
      { // фаза 3: сброс скорости без рывка до a*a/(2*Rd) c a=-Ad
        r = 0.;
        //!!!if (2. * (v + a) * R < a2)
        //!!!  r = 2. * (a2 / (2. * R) - v - a);
        RFILT_DBG("stop phase 3");
      }
      else // a > -Ad && 2*v*Rd > a2 && a<=0 && v<=v+(a2-Ad2)/(2*Rd)
      { // фаза 2: набор/добор тормозного ускорения до -Ad
        r = -R;
        if (a + r < -self->Ad)
          r = -self->Ad - a;
        RFILT_DBG("stop phase 2");
      }
    }
#else
    r = -R;
    RFILT_DBG("stop phase");
#endif
  }
  else // if (s >= 0. && s <= dx)
  { // 0 <= s < dx - движение в нужном направлении
    // v>0 a<0 - торможение (как это?)
    // v<0 a>0 - выход из реверса
    // v<0 a<0 - провал в реврс
    // v>0 a>0 - разгон
    if (v < 0.)
    { // реверс
      r = R;
      RFILT_DBG("reverse v<0");
    }
    else // v >= 0
    { // v >= 0, 0 <= s <= dx
#if 0 // FIXME: эффективность этого кода не велика, но эффект есть!
      // тормозить "рано", но проверить не "поздно" ли будет на следующем такте
      double s1, v1, a1, nt;
      double r1 = rfilt_limit_r(self, v, a, +R); // рывок разгона
      r         = rfilt_limit_r(self, v, a, -R); // рывок торможения
      s1 = v + a * 0.5 + r1 / 6.;
      v1 = v + a + r1 * 0.5;
      a1 = a + r1;
      s1 += rfilt_s(self, R, v1, a1, &nt);
      if (s1 >= dx)
      { // на следующем такте тормозить уже поздно =>
        // пропорциональное управление рывком на текущем такте:
        //   r =limit(-R) -> путь s
        //   r1=limit(+R) -> путь s1
        if (s1 > s)
          r = (dx - s) / (s1 - s) * (r1 - r) + r;
        else
          r = 0.;
      }
      else
        r = R;
#else
      r = R;
#endif

#if 0 // FIXME: вероятность выполнения этого кода ничтожна
      if (dx == 0. && a == 0.)
      { // частный (вырожденный) случай s=dx=0, a=0, v=0
        r = 0.; // никуда не разгоняемся - покой
      }
#endif
      RFILT_DBG("accelerate");
    }
  }

  // ограничить рывок с учетом ограничения скорости и ускорения
  r = rfilt_limit_r(self, v, a, r);

  // восстановить знаки
  return r * sig;
}
#endif
//-----------------------------------------------------------------------------
// расчёт рывка
double rfilt_r(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double dx)     // невязка заданой и текущей позиции
{
  double r1, r2, s1, s2, v1, v2, a1, a2, r, nt1, nt2;
  double R = self->R;
  double sig = (dx >= 0.) ? 1. : -1.;

  // привсти знаки относитльн dx
  v  *= sig;
  a  *= sig;
  dx *= sig; // теперь dx >= 0
  
  dx *= (1. - self->B); // FIXME

  r1 = rfilt_limit_r(self, v, a, -R);
  r2 = rfilt_limit_r(self, v, a,  R);
  s1 = v + a * 0.5 + r1 / 6.;
  s2 = v + a * 0.5 + r2 / 6.;
  v1 = v + a + r1 * 0.5;
  v2 = v + a + r2 * 0.5;
  a1 = a + r1;
  a2 = a + r2;
  s1 += rfilt_s(self, self->R, v1, a1, &nt1);
  s2 += rfilt_s(self, self->R, v2, a2, &nt2);

  if (dx <= s1 && dx <= s2)
  { // торможение
    double a2 = a * a;
    double Ad2 = self->Ad * self->Ad;
    if (s1 <= s2) { r = r1; self->s = s1; self->nt = nt1; } // dx <= s1 <= s2 
    else          { r = r2; self->s = s2; self->nt = nt2; } // dx <= s2 <  s1
    // торможение состоит из от 1-й до 4-х фаз:
    // 1. сброс ускорения разгона с рывком r<0 и набором сокрости (если a>0)
    // 2. набор (добор) тормозного ускорения со сбросом скорости и с рывком r<0
    // 3. сброс скорости с постоянным ускорением без рывка r=0 (если ускорение достигло максимума)
    // 4. сброс скорости и тормозного ускорения до нуля с положительным рывком r>0
    if (a > 0.)
    { // фаза 1: требуется сброс ускорения разгона до 0 и далее
      r = -R;
      RFILT_DBG("stop phase 1");
    }
    else // if a <= 0
    {
      if (2. * v * R <= a2)
      { // фаза 4: требуется сброс тормозного ускорения и остаточной скорости до нуля
        //r = R;
        r = a2 / (2. * v); // FIXE
        if (r > -a) r = a;
        // FIXME:
        //!!!self->v = a2 / (2. * R); // коррекция накапливаемых ошибок вычислений!
        RFILT_DBG("stop phase 4");
      }
      else if (v <= v + (a2 - Ad2) / (2. * R) ||
               a <= -self->Ad)
      { // фаза 3: сброс скорости без рывка до a*a/(2*Rd) c a=-Ad
        r = 0.;
        //!!!if (2. * (v + a) * R < a2)
        //!!!  r = 2. * (a2 / (2. * R) - v - a);
        RFILT_DBG("stop phase 3");
      }
      else // a > -Ad && 2*v*Rd > a2 && a<=0 && v<=v+(a2-Ad2)/(2*Rd)
      { // фаза 2: набор/добор тормозного ускорения до -Ad
        r = -R;
        if (a + r < -self->Ad)
          r = -self->Ad - a;
        RFILT_DBG("stop phase 2");
      }
    }
    r = rfilt_limit_r(self, v, a, r);
  }
  else if (dx >= s1 && dx >= s2)
  { // разгон
    if (s1 >= s2) { r = r1; self->s = s1; self->nt = nt1; } // s2 <= s1 <= dx
    else          { r = r2; self->s = s2; self->nt = nt2; } // s1 <  s2 <= dx
    RFILT_DBG("accelerate (Rmax)");
  }
  else
  { // пропорциональное управление рывком
    if (s1 <= s2) { self->s = s1; self->nt = nt1; } // s1 < dx < s2
    else          { self->s = s2; self->nt = nt2; } // s2 < dx < s1
    
    if (s1 == s1)
      r = 0.;
    else
      //r = (dx * (1. - self->B)- s1) * (r2 - r1) / (s2 - s1) - r1;
      r = (dx - s1) * (r2 - r1) / (s2 - s1) - r1;
    RFILT_DBG("accelerate (calc R)");
  }
  
  // восстановить знаки
  self->s *= sig;
  return r * sig;
}
//-----------------------------------------------------------------------------
// проверка возможности перемещения в заданноую позицию за один такт
// (если 1, то можно осуществить прямой переход и обнулять текущие a и v)
static inline int rfilt_jump(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double dx)     // невязка заданой и текущей позиции
{
#if 0
  if (fabs(dx) > self->R / 6.)
    // превышен предельно максимальный путь до покоя за такт
    // даже без учета a и v
    return 0;
#endif

#if 1
  if (fabs(dx) > self->R / 24.)
    // превышен максимальный путь из покоя в покой за такт
    // даже без учета a и v
    return 0;
#endif
  
  return 1;
}
//-----------------------------------------------------------------------------
#if 0 // FIXME: эксперментальный код
// проверка возможности перемещения в заданноую позицию за один такт
// (если 1, то можно осуществить прямой переход и обнулять текущие a и v)
int rfilt_jump_exp(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра (константы)
  double v,      // приведенная текущая скорость
  double a,      // приведенное текущее ускорение
  double dx)     // невязка заданой и текущей позиции
{
  double sa, sb, na, nb, m;
  double sig = (dx >= 0.) ? 1. : -1.;

  // вычислить тормозной путь и его время
  sa = rfilt_s(self, self->R, v, a, &na);

  // вычислить обратный тормозной путь и его время
  sb = rfilt_s(self, self->R, -v, a, &nb);

  // привсти знаки относительно dx
  v  *= sig;
  a  *= sig;
  sa *= sig;
  sb *= sig;
  dx *= sig; // теперь dx >= 0

  if (sa > dx)
    return 0; // перебор тормозного пути
    
  if (na > 1)
    return 0; // даже затормозить не успеваем
    
#if 0
  if (fabs(dx) > self->R / 6.)
    return 0; // превышен предельно максимальный путь в покой за такт
#endif

#if 1
  // FIXME: следующяя проверка приблизительная!
  m = pow(24. * fabs(dx + sb) / self->R, 1./3.) - nb;
  if (m > 1.)
    return 0; // превышена величина масимального пути из покоя за один такт
#endif
  
  return 1;
}
#endif
//-----------------------------------------------------------------------------
// выполнить шаг фильтрации
// (функция возвращает "профильтрованное" значение с учетом ограничений)
double rfilt_step(
  rfilt_t *self, // указатель на внутреннюю структуру фильтра
  double y)      // входное значение для фильтра
{
  double dx, r, xs;

  // ограничить заданное значение
  y = RFILT_LIMIT(y, self->Xmin, self->Xmax);

  // требуемое перемещение
  dx = y - self->x;

#if 0 // FIXME
  // проверить возможность "малого хода"
  if (rfilt_jump(self, self->v, self->a, dx))
  {
    RFILT_DBG("pipe mode: x=y=%f dx=%f", y, dx);
    self->v = self->a = 0.;
    self->s = self->nt = 0.;
    return self->x = y;
  }
#endif

  // вычислить рывок
  r = rfilt_r(self, self->v, self->a, dx);

  xs = self->x + self->s;

  RFILT_DBG("work: y=%f x=%f v=%f a=%f r=%f dx=%f s=%f x+s=%f nt=%f",
             y, self->x, self->v, self->a, r, dx, self->s, xs, self->nt);

  // выполнить расчёт следующей позиции, скорости и ускорения
  self->x += self->v + self->a * 0.5 + r / 6.;
  self->v += self->a + r * 0.5;
  self->a += r;

  return self->x;
}
//-----------------------------------------------------------------------------

/*** end of "rfilt.c" file ***/
