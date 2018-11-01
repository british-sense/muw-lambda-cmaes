# (*μ*<sub>W</sub>, *λ*)-CMA-ES Algorithm


C++における(*μ*<sub>W</sub>, *λ*)-CMA-ESの実装です.
実装したCMA-ESは以下の論文のものです.


Completely Derandomized Self-Adaption in Evolution Strategies
http://www.cmap.polytechnique.fr/~nikolaus.hansen/cmaartic.pdf


以下の環境で実験をしました.

|parameter|value|
|:-:|:-:|
|λ|10|
|μ|2|
|w_i|1|
|σ|100|
|n|100|
|generation|1000|
|benchmark|rastrigin|

実験結果

![rastrigin1](https://github.com/ko-cha/CMA-ES/blob/master/img/image002.png "rastrigin1")

![rastrigin2](https://github.com/ko-cha/CMA-ES/blob/master/img/image004.png "rastrigin2")

100次元でのrastriginにおいて最適解を得ることができました.
