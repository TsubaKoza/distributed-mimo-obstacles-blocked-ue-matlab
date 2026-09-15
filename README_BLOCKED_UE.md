# 単独遮蔽UEの協調送信・電力全探索

## 実行

MATLABでこのコードフォルダをcurrent folderにする。必須製品はbase MATLABのみ。

```matlab
runRayTracingSanityChecks
runBlockedUEPowerSanityChecks
experiment = main_blocked_ue_power_experiment();
```

```matlab
cfg = configBlockedUEExperiment(struct( ...
    'M',20,'K',80,'N_AP',4,'powerStep',0.05, ...
    'maxCandidateAPsForPowerSweep',4,'targetBlockedUEIndex',[], ...
    'showFigures',false,'verifyReproducibility',true));
experiment = main_blocked_ue_power_experiment(cfg);
```

`targetBlockedUEIndex` は元のUE番号（遮蔽UEリスト内の順位ではない）。指定UEが全条件を満たさなければ採用しない。各attemptのseedと採用条件をCSVに残し、最大試行数で停止する。途中でLoSを削除したり、UE周囲に建物を追加したりしない。

既存処理は `main_ray_tracing_3D()` で従来通り実行できる。`random3D` のAP/UE/建物生成順と乱数の消費順を維持。従来のランダム複素反射係数とpruningなしが、この入口のデフォルトである。都市配置だけをRTする場合も `main_ray_tracing_3D(configBlockedUEExperiment())` が利用できる。

## 初期条件

|設定|デフォルト|
|---|---|
|搬送周波数 `fc`|100e9 Hz = 100 GHz|
|波長・素子間隔|3 mm・1.5 mm|
|領域|500 m × 500 m、高さ表示45 m|
|AP / UE / AP当たりアンテナ|`M=20` / `K=80` / `N_AP=4`|
|AP / UE高さ|8 m / 1.5 m (`zAP`, `zUE`)|
|AP最小距離|50 m (`minimumAPSeparation`)|
|縦横corridor数|各3本 (`urbanCorridorCount`)|
|corridor幅|18 m (`urbanCorridorWidth`)|
|建物幅 / 奥行 / 高さ候補|20–60 / 20–60 / 10–40 m|
|建物間の最小空間|8 m (`urbanBuildingGap`)|
|yaw jitter|0度、設定可能（最大15度）|
|建物沿いのUE比率|0.65 (`urbanNLoSProneUEFraction`)|
|反射係数|固定0.6∠0度|
|pruning|NLoS scalar path power ≥ −160 dB、最大本数無制限|
|serviceable blocked UE必要数|5 (`minBlockedUECount`)|
|最大シナリオ試行|10 (`maxScenarioGenerationAttempts`)|
|対象UEの使用可能NLoS AP数|3以上|
|帯域幅 / 目標速度|20 MHz / 10 Mbps|
|thermal noise density / NF|−174 dBm/Hz / 7 dB（暫定値）|
|AP当たり上限 / 刻み|1 W / 0.05 W|
|sweep対象AP|実効利得上位4台まで|
|最大組合せ数 / 保存上位解|1,000,000 / 50|
|再現性の再計算|有効（geometryと全リンクRTをもう一度計算）|

都市は縦横corridor間の16街区をさらに矩形区画へ分割する。デフォルトでは64棟。幅・奥行は指定範囲から採り、回転後の外接矩形が区画に収まるものを採用するため、実際の寸法分布は区画サイズで条件付けられる。領域やcorridor数を変更すると区画数・棟数も変化する。都市の建物数は`numObstacles`では指定せず、街区と寸法から決まる。`numObstacles` は従来random3D用の意味を保持する。

各区画間には通り抜けられる隙間があり、UEを囲うための壁は置かない。APはcorridorの両側に分散配置。UEは街路と全方向の建物外周からサンプリングし、APに対する可視性を調べて配置を移動しない。corridor、建物区画、建物沿いという配置ラベルはLoSフラグを決めない。

## 既存チャネルとMRTの厳密な対応

既存READMEの定義に従い、保存列ベクトルを

\[
g_m=\texttt{result.h\_true(:,m,k)}\in\mathbb C^{N_{AP}}
\]

と書く。これはUE→APの受信チャネルであり、既存の計算は

\[
g_m=\sum_{\ell\in\mathcal P_{mk}}\alpha_\ell a_\ell,
\quad a_{\ell,n}=e^{-j2\pi\Delta r_n^T u_\ell/\lambda}.
\]

`alpha=complexGain` は伝搬位相を含み、LoS/反射ではFriis field amplitude、回折では既存knife-edge損失を用いる。新しいMRTでpathごとに位相を揃え直さない。path同士の強め合い・打ち消し合いはそのまま残る。

今回明示するdownlink前提は**相反な狭帯域チャネル、正確な複素チャネル情報、AP間の完全な位相・時間同期**である。UE→APベクトルを使ったdownlinkは `g.'*w`（非共役転置）。通常のHermitian表記との対応のためだけに \(h_m^{DL}=g_m^*\) と置くと、

\[
(h_m^{DL})^H w_m=g_m^T w_m,
\quad w_m=\frac{g_m^*}{\|g_m\|_2},
\quad \|w_m\|_2^2=1.
\]

```matlab
g = result.h_true(:,m,k);
w = conj(g)/norm(g);
hDL = conj(g);
a = g.'*w;        % same as hDL'*w; do NOT use g'*w here
G = abs(a)^2;
```

ゼロベクトルではweightと実効利得を0とし、電力候補から外す。

\[
a_m=g_m^T w_m=\|g_m\|_2,\qquad G_m=|a_m|^2=\|g_m\|_2^2.
\]

従ってこのunit-norm MRT baselineではChannelGainとEffectiveGainは一致する。チャネル利得は既存の `sum(abs(h_true).^2)` であり、素子平均やpath電力和へ置き換えない。ランキングのStrongestPathPowerは `10*log10(abs(complexGain)^2)`（単一pathのscalar field gainの電力、送信電力は未適用）。配列合成利得とは区別する。

## coherent joint transmission・雑音・速度

各APは同じシンボル \(s\)、\(E|s|^2=1\) を送る。AP単位の電力を \(p_m\) として

\[
x_m=\sqrt{p_m}w_m s,\quad E\|x_m\|^2=p_m,
\quad y=\left(\sum_m\sqrt{p_m}a_m\right)s+n,
\]

\[
P_{desired}=\left|\sum_m\sqrt{p_m}a_m\right|^2,
\quad P_{total}=\sum_m p_m,\quad 0\le p_m\le P_{max}.
\]

**電力をreflection ray / diffraction rayへ割り当てない。** AP間の電力和ではなく、複素振幅を合計した後で絶対値二乗を取る。AP間位相ずれ、hardware impairment、帯域内の周波数変動はこのbaselineに含めない。

\[
N_{dBm}=-174+10\log_{10}(B)+NF_{dB},
\quad \sigma^2=10^{(N_{dBm}-30)/10}\;[W].
\]

20 MHz、NF=7 dBでは**−93.9897 dBm、3.99052463×10⁻¹³ W**。NFは比較用の暫定値であり、100 GHz受信機に普遍的な値ではない。100 GHzの搬送周波数を雑音帯域幅として使わない。

**single blocked UE experimentのため他UE干渉は存在せず、今回のSINRは実質SNRである。**

\[
\mathrm{SINR}=\mathrm{SNR}=P_{desired}/\sigma^2,
\quad R=B\log_2(1+\mathrm{SINR}),
\quad \gamma=2^{R_{target}/B}-1.
\]

デフォルトでは \(\gamma=\sqrt2-1\simeq0.414214\)、−3.82776 dB。毎回設定から再計算する。これは狭帯域Shannon baselineの速度であり、符号化・制御オーバーヘッド込みの実効スループットではない。

## pruning・遮蔽判定・対象選択

元コードには `pruneRayPaths.m` はなかった。追加したpruningはNLoS pathのscalar gain閾値と本数上限であり、**幾何的LoSを削除しない**。旧入口では無効なので元のチャネルを変えない。新実験では閾値−160 dB。本閾値は研究用の数値カットであり、受信機感度・QoS判定ではない。感度は送信電力と雑音にも依存する。

全APの `hasLoS=false`、pruning後のreflection/diffractionがあるAPが2台以上、全AP outageでないことをserviceable blocked UE条件とする。`outageMask` は既存通り「pathが0本」の意味で、QoS未達の意味へ変更しない。

対象は、さらに非ゼロ合成利得のNLoS APが3台以上で、sweepに使う上位AP各1 Wで設定QoSを満たすもの。自動選択はNLoS AP数の降順、同数なら最大sweep速度の降順、さらに同数ならUE番号順。明示指定でも全条件を検査する。シナリオ拒否は遮蔽UE数不足または適切な対象不在で行い、最大試行数に達すると診断MAT/CSVを保存して停止する。

## 全探索・最小電力・相関

0から上限までの格子（端点上限を必ず含む）を使う。0.05 W刻み・1 W上限・4 APなら21⁴=194,481通り。1 Wに割り切れない刻みの場合は最後の区間だけ短くなる。組合せ数を事前表示し、上限超過時はwarningの後に停止する。全使用可能APのランキングを保存するが、最適解は**選択した上位AP・離散格子内の最小値**。全APや連続電力上の大域最適と呼ばない。

単独APの連続必要電力は

\[
p_{m,min}=\gamma\sigma^2/G_m.
\]

1 W超のAPも表・相関に残す。10 Mbpsの固定列は目標速度設定を変えても10 Mbpsの意味を保ち、設定目標に対応する別列 `MinimumPowerForTargetRate_W` を併記する。`FeasibleWithin1W` と `AchievedRateAt1W_Mbps` も文字通り1 W、変更可能なAP上限用には別列を使う。

QoS達成表は合計電力の昇順、同電力なら速度の降順、同速度ならCombinationID順。上位50解と最良解のAP別電力・配分率を保存する。全ゼロ電力も全探索に含まれる。

Pearsonは実効利得[dB]と単独10 Mbps必要電力[W]、Spearmanは同じ量の平均順位を使用。有限な標本数Nも出力する。Pearsonのp値はbase MATLABの`betainc`で両側t参照分布を計算する。SpearmanはN≤8なら全ラベル順列による両側の正確な置換p値、N>8ならt近似。N<3や定数系列ではp値をNaNにする。APリンクは選択・幾何的依存があるため通常の独立標本仮定は保証されない。

**必要電力はモデル上利得の逆数なので、Spearmanの−1は通常数学的に決まる。これだけで未知の物理現象や統計的因果を実証したとはいえない。** 小標本・シナリオのrejection sampling・対象選択による条件付けもあり、この1 UEから都市全体へ一般化しない。

## 保存物

`blocked_ue_results/blocked_ue_experiment.mat` は `scenario`, `targetBlockedUE`, `channelRankingTable`, `singleAPMinimumPowerTable`, `allPowerCombinationTable`, `feasiblePowerCombinationTable`, `minimumPowerSolution`, `correlationResults`, `config` をトップレベル変数として保存。加えて全RT結果、分類、MRT weights、実効振幅、上位解、対象候補表、試行履歴、requestedConfig、validationを保存する。

各tableのCSV、最良解CSV、相関CSV、実測値をまとめた `RUN_SUMMARY.md` と5種類のPNG/FIGも生成する。

```matlab
e = load('blocked_ue_results/blocked_ue_experiment.mat');
e.channelRankingTable
e.singleAPMinimumPowerTable
e.minimumPowerSolution.APAllocation
e.topPowerCombinationTable
validateBlockedUEExperiment(e,true)  % 同seedの全RTも再計算
```

## 物理モデルの限界・反射係数

**今回の単独遮蔽UE実験では、遮蔽物材質差による影響を排除するため、全遮蔽物で同じ複素反射係数を使用している。** `reflectionCoefficientMode="fixed"`, `fixedReflectionMagnitude=0.6`, `fixedReflectionPhaseDeg=0`。全障害物で同値をassertする。比較用に`"random"`を選べるが都市研究のデフォルトは必ずfixed。

既存Ray Tracingのarray response、伝搬位相と回折近似はそのまま用いた。保存チャネルを新しい受信方程式へ接続するため、前述の相反性・完全CSI・位相同期を明示した。実測アンテナやRFチェーンの校正モデルは存在せず、今回は追加しない。既存array responseの符号をここで再解釈して反転することも行わない。

特に、既存回折は障害物を局所的に零厚knife edgeへ置換し、当該障害物自体の交差を無視して2区間の見通しを検査する。厚い建物の厳密UTDではなく、幾何的に建物内部を横切るように見える折れ線を許し得る。ITU型式は振幅損失のみで追加回折位相を含まず、折れ線長による位相だけを使う。従って協調送信の数値は**現行の近似モデル内での検証**であり、そのまま実都市の必要電力予測として主張しない。

材質・偏波・入射角依存のFresnel係数、大気吸収、厳密UTD、二回以上の反射/回折、実機同期誤差、beam squint、周波数選択性は未実装。今回は既存物理モデルを維持する要件のため追加しない。複数UE同時通信・他UE干渉・GE-MCTS・推定・pilot・OFDM等も追加していない。

## 変更一覧と検証

変更した既存ファイルは `configRayTracing3D.m`（追加設定）、`generateScenario3D.m`（都市分岐）、`generateObstacles3D.m`（係数モード）、`generateRayTracingChannels3D.m`（pruningと残存数）、`README.md`（案内）。幾何的LoS・反射・回折の探索関数、AP配列、既存main、既存sanityは変更していない。

|追加ファイル|役割|
|---|---|
|configBlockedUEExperiment.m|研究用の初期値・入力検証|
|main_blocked_ue_power_experiment.m|再生成・対象選択・計算・保存の入口|
|generateUrbanStreetCanyon3D.m|街区と屋外AP/UE生成|
|validateUrbanScenario3D.m|重なり・corridor・屋外・間隔・固定係数検査|
|sampleReflectionCoefficient.m|fixed/random係数（旧乱数順維持）|
|pruneRayPaths.m|任意のNLoS経路pruning|
|identifyServiceableBlockedUEs.m|全AP遮蔽とserviceable判定|
|selectTargetBlockedUE.m|対象UE条件・選択理由|
|computeBlockedUEMRT.m|正規化MRTと実効チャネル|
|computeExperimentNoise.m|熱雑音・NF・必要SINR|
|buildBlockedUEChannelRanking.m|全APランキング|
|evaluateCoherentPowerAllocations.m|複素振幅合成・SINR・速度|
|sweepBlockedUEPower.m|全組合せ・単独必要電力・最良/上位解|
|analyzeGainPowerCorrelation.m|相関・p値・標本数|
|plotBlockedUEExperiment.m|5種類の図|
|saveBlockedUEExperiment.m|MAT/CSV/実行サマリー保存|
|validateBlockedUEExperiment.m|実験結果のassert・同seed再現性|
|runBlockedUEPowerSanityChecks.m|解析式・特殊ケース・旧sanity・旧結果回帰|

同seed回帰は変更前のMAT snapshotを渡すことで再検査できる。

```matlab
runBlockedUEPowerSanityChecks('path/to/baseline.mat')
```

このsnapshotは現在のコードから生成し直したものではなく、必ず変更前コードで作成したものを使用する。


## 電力刻みを評価するための連続電力参考解

`computeContinuousPowerReference.m` は、同じsweep候補APに対する連続電力の最小値を別途計算する。総当たり解を置き換えない。`continuousPowerReference` とCSVに保存する。

`u_m=sqrt(p_m)`、`a_m=sqrt(G_m)` と置くと、目的は `sum(u_m^2)` の最小化、制約は `sum(a_m*u_m)>=sqrt(gamma*noise)`、`0<=u_m<=sqrt(Pmax)`。これは凸二次最小化であり、KKT条件から `u_m=min(sqrt(Pmax),lambda*a_m)`。1変数の単調な制約を二分法で解く。

AP上限に達しない場合は `p_m=gamma*noise*G_m/(sum(G))^2`、`Ptotal=gamma*noise/sum(G)`。従って、この理想的同期モデルでAP回路電力等を考えなければ、強いリンクほど大きい電力を与えつつ複数APに分散することが連続電力では有利になる。0.05 W刻みの最適解がAP単独利用だったとしても、連続電力でも単独利用が最適という結論にはならない。
