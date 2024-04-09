simulation
================
2024-04-09

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.1      ✔ purrr   0.3.4 
    ## ✔ tibble  3.2.1      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.0      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.2      ✔ forcats 0.5.2

    ## Warning: 패키지 'ggplot2'는 R 버전 4.2.3에서 작성되었습니다

    ## Warning: 패키지 'tibble'는 R 버전 4.2.3에서 작성되었습니다

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(MASS)
```

    ## 
    ## 다음의 패키지를 부착합니다: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(tictoc)
library(GIGrvg)
```

    ## Warning: 패키지 'GIGrvg'는 R 버전 4.2.3에서 작성되었습니다
