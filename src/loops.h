#ifndef LOOPS_H
#define LOOPS_H

#define LOOP_KJI(LOOP_DA,LOOP_EXPRESSION) \
        do {                                                            \
                int LOOP_XS,LOOP_YS,LOOP_ZS,LOOP_XM,LOOP_YM,LOOP_ZM;    \
                DMDAGetCorners(LOOP_DA,&LOOP_XS,&LOOP_YS,        \
                                      &LOOP_ZS,&LOOP_XM,&LOOP_YM,       \
                                      &LOOP_ZM);          \
                for (int k=LOOP_ZS; k<LOOP_ZS+LOOP_ZM; k++) {           \
                        for (int j=LOOP_YS; j<LOOP_YS+LOOP_YM; j++) {   \
                                for (int i=LOOP_XS; i<LOOP_XS+LOOP_XM;  \
                                     i++) {LOOP_EXPRESSION;}}}} while(0)

#endif /* LOOPS_H */
