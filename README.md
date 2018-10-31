# CMA-ES


C++におけるCMA-ESの実装です.
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
![rastrigin1](https://github.com/ko-cha/CMA-ES/img/image002.png "rastrigin1")
![rastrigin2](https://github.com/ko-cha/CMA-ES/img/image004.png "rastrigin2")
