# ハッシュ関数で利用される同種写像グラフの効率的な連結性判定法
このリポジトリは卒論のソースコードを含みます. 

`run_experiments.py `
で以下の値を設定する

### 設定項目
|  パラメータ | 設定値の例 | 説明 |
| :--- | :--- | :--- |
| `MODE` | `"ssp"` / `"ssg"` | 同種写像グラフの種類 (超特別 / 超特異非超特別) を選択 |
| `PRIME_MIN` | $7$ | 実行する素数の最小値 |
| `PRIME_MAX` | $100$ | 実行する素数の最大値 |
| `MAX_WORKER` | $1$ | 並列実行するプロセス数 |

```
python3 run_experiments.py 
```
を実行すると, この例では, について超特別 / 超特異非超特別の連結性を素数 $p\in[7,100]$ について $1$並列で判定する.

各素数

`data/ssp_results`または`data/ssg_results`に保存される
| 項目 | 意味 | 
| :--- | :--- | 
| p | 基礎体の標数 |
| d/d1 |  $\mathbb{F}_p$ のある平方非剰余な元, 二次拡大$\mathbb{F}_{p^2} := \mathbb{F}_p[x] / (x - d^2)$ |に用いる
|d1c0, d2C1|
|MODE |supr||



p,d1,d2c0,d2c1,construction_time_sec,jacobian_count,vertex_count,edge_count,connectivity_check_time_sec,is_strongly_connected,scc_count,are_all_jacobians_connected_via_main_scc,are_all_jacobian_pairs_connected_via_good_paths,are_all_minor_sccs_size_1,are_all_minor_nodes_self_isogenies,all_minor_sccs_have_preds_from_main,any_minor_sccs_have_succs_to_main,exists_minor_to_minor_edge,condensed_graph,construction_peak_bytes,connectivity_check_peak_bytes