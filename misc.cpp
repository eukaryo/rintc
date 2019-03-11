/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"misc.h"
#include"parameter.h"

namespace rintdwr {

int Ceiling2Power(int h) {

	//Returns the minimum integer as represented by the power of 2,
	//greater than or equal to h.

	//h以上の整数のうち、2のべき乗で表されるような最小の数を求めて返す。

	assert(0 <= h && h <= 1000000);

	--h;
	h |= h >> 1;
	h |= h >> 2;
	h |= h >> 4;
	h |= h >> 8;
	h |= h >> 16;
	return h + 1;
}

std::string MakeRandomRnaSeq(const int length, const int seed) {
	std::mt19937_64 rnd(seed);
	std::uniform_int_distribution<int> kind(0, 3);
	const std::string base("ACGU");
	std::string sequence;
	for (int i = 0; i < length; ++i)sequence.append(1, base[kind(rnd)]);
	return sequence;
}

int ComputeMaxLoop(const std::string& structure) {
	//Compute the maximum value of unpaired base number
	//of "internal loop" and "bulge" in "structure".

	//RNA二次構造structureのなかで、internal loopとbulgeのunpaired base数の最大値を求めて返す。

	const int n = int(structure.size());

	std::vector<int> bp(n, -1);
	std::stack<int> stk;
	for (int i = 0; i < n; ++i) {
		switch (structure[i]) {
		case '(':
			stk.push(i);
			break;
		case ')':
			bp[i] = stk.top();
			bp[stk.top()] = i;
			stk.pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}

	int ans = 0;
	for (int pos = 1; pos < n - 1; pos++) {
		if (structure[pos] != '.')continue;
		int dangleflag = 1;
		int danglenum = 1;
		int branchnum = 0;
		for (int i = pos + 1; i != pos; ++i) {
			if (bp[i] != -1) {
				++branchnum;
				i = bp[i];
			}
			else ++danglenum;
			if (i == n - 1) {
				danglenum = 0;
				break;
			}
		}
		if (branchnum == 2)ans = std::max(ans, danglenum);
	}
	return ans;
}

std::vector<std::vector<int>> VerifyAndParseStructure(const std::string& structure, const std::string& sequence, const int max_span, const int max_loop) {
	//引数structureはブラケット表記のRNA二次構造とする。
	//これのベリファイをしつつ、何番と何番の塩基が塩基対を組んでいるか調べて返す。

	//返り値は行列で、i<jなる(i,j)が塩基対を組んでいるなら[i][j]=1で、それ以外は0とする。
	//このときのi,jは1-originとする。

	for (const char c : sequence) {
		assert(c == 'A' || c == 'U' || c == 'G' || c == 'C');
	}
	assert(structure.size() == sequence.size());
	assert(1 <= max_span && max_span <= int(structure.size()));

	const int n = int(sequence.size());
	const std::string bp = "AU UA GC CG GU UG";
	std::string query = "XX";
	std::vector<std::vector<int>>ans(n + 1, std::vector<int>(n + 1, 0));
	std::stack<int> bp_pos;
	for (int i = 1; i <= n; ++i) {
		switch (structure[i - 1]) {
		case '(':
			bp_pos.push(i);
			break;
		case ')':
			assert(bp_pos.size() >= 1);
			assert(TURN < (i - bp_pos.top()) && (i - bp_pos.top()) <= max_span);

			query[0] = sequence[bp_pos.top() - 1];
			query[1] = sequence[i - 1];
			assert(bp.find(query) != std::string::npos);

			ans[bp_pos.top()][i] = 1;
			bp_pos.pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}
	assert(bp_pos.size() == 0);
	assert(ComputeMaxLoop(structure) <= max_loop);
	return ans;
}

std::vector<std::vector<int>> ComputePredistanceMatrix(const std::vector<std::vector<int>>& S) {
	//引数SはVerifyAndParseStructureで得られた塩基対行列S(1-origin)で、
	//sequenceの二次構造としてvalidだとする。

	//PredistanceMatrix(1-origin)を作って返す。
	//PredistanceMatrixとは[Mori et al., 2014]のsuppの式(S10)のCのことである。

	const int n = int(S.size()) - 1;
	std::vector<std::vector<int>>C(n + 1, std::vector<int>(n + 1, 0));

	//[Mori et al., 2014]のsuppのsection S2ではDPで記述されているが、
	//ここではメモ化再帰で書く。
	std::vector<std::vector<int>>memo(n + 1, std::vector<int>(n + 1, 0));

	//[Mori et al., 2014]のsuppの式(S11)、(S12)の初期化処理である。
	for (int i = 1; i <= n; ++i) {
		memo[i][i] = 1;
		C[i][i] = 0;
	}
	for (int i = 1; i <= n - 1; ++i) {
		memo[i][i + 1] = 1;
		C[i][i + 1] = S[i][i + 1];
	}

	//[Mori et al., 2014]のsuppの式(S13)の再帰処理である。
	std::function<int(int, int)> GetC = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (memo[i][j])return C[i][j];
		memo[i][j] = 1;
		return C[i][j] = S[i][j] + GetC(i + 1, j) + GetC(i, j - 1) - GetC(i + 1, j - 1);
	};

	//メモ化再帰による計算を実行する。
	for (int i = 1; i <= n; ++i) {
		for (int j = i; j <= n; ++j) {
			GetC(i, j);
		}
	}

	return C;
}

int ComputeMaxHammingDistance(const std::string& sequence, const std::vector<std::vector<int>>& S, const int max_span, const int max_loop) {
	//引数sequenceはAUCGから成るRNA配列とする。
	//引数SはVerifyAndParseStructureで得られた塩基対行列S(1-origin)とする。

	//sequenceについて可能な全ての二次構造についての
	//reference_structureからのハミング距離の最大値を求めて返す。
	//これは[Mori et al., 2014]のsuppのsection S3の処理である。

	//ハミング距離の最大値を求めるだけならNassinov型DPでもできるが、ここではあえてMcCaskill型DPを使っている。理由は、
	//本番のMcCaskill型DP計算で対象とされる二次構造集合と、完全に同じ集合(maxloop制約を含む)を扱いたいからである。

	const int n = int(sequence.size());
	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(S);

	//メモ化再帰のためのデータ構造である。firstは初期化時はfalseで、値を計算したらtrueにして、secondに値を入れる。
	typedef std::pair<bool, int>MemInt;

	//[Mori et al., 2014]のsuppの式(S23)～(S27)で定義されているDである。
	//std::vectorの初期化が0埋めでなされるので式(S23)がなされたとみなす。
	std::vector<std::vector<MemInt>> D(n + 1, std::vector<MemInt>(n + 1, std::make_pair(false, 0)));
	std::vector<std::vector<MemInt>> D1(n + 1, std::vector<MemInt>(n + 1, std::make_pair(false, 0)));
	std::vector<std::vector<MemInt>> Db(n + 1, std::vector<MemInt>(n + 1, std::make_pair(false, 0)));
	std::vector<std::vector<MemInt>> Dm(n + 1, std::vector<MemInt>(n + 1, std::make_pair(false, 0)));

	//[Mori et al., 2014]の式(19)～(23)の左辺のZたちである。
	//[Mori et al., 2014]のnotationに合わせて添字を1-originとする。
	std::function<int(int, int)> GetD;
	std::function<int(int, int)> GetD1;
	std::function<int(int, int)> GetDb;
	std::function<int(int, int)> GetDm;

	GetD = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return 0;
		if (j == i - 1)return 0;
		if (D[i][j].first)return D[i][j].second;

		int m = 0;
		for (int k = i - 1; k < j; ++k) {
			m = std::max(m, GetD(i, k) + GetD1(k + 1, j) - C[i][k] - C[k + 1][j]);
		}

		return (D[i][j] = std::make_pair(true, m + C[i][j])).second;
	};

	GetD1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return -10000000;//[Mori et al., 2014]のSuppの式(S23)ではExactな値が求まらない(余計に大きい値が求まることがある)
		if (D1[i][j].first)return D1[i][j].second;

		int m = -10000000;
		for (int k = i; k <= j && (k - i) <= max_span; ++k) {
			m = std::max(m, GetDb(i, k) - C[i][k]);
		}
		return (D1[i][j] = std::make_pair(true, m + C[i][j])).second;
	};

	GetDb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return -10000000;//[Mori et al., 2014]のSuppの式(S23)ではExactな値が求まらない(余計に大きい値が求まることがある)
		if (Db[i][j].first)return Db[i][j].second;

		//iとjが塩基対を組み得ないならゼロを返す。
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || !(i + TURN < j) || (j - i) > max_span)return (Db[i][j] = std::make_pair(true, -10000000)).second;

		int m = 0;
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				m = std::max(m, GetDb(k, l) - C[k][l]);
			}
		}
		for (int k = i + 2; k <= j - 1; ++k) {
			m = std::max(m, GetDm(i + 1, k - 1) + GetD1(k, j - 1) - C[i + 1][k - 1] - C[k][j - 1]);
		}
		return (Db[i][j] = std::make_pair(true, m + C[i][j] - 2 * S[i][j] + 1)).second;
	};

	GetDm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return -10000000;//[Mori et al., 2014]のSuppの式(S23)ではExactな値が求まらない(余計に大きい値が求まることがある)
		if (j == i - 1)return -10000000;//[Mori et al., 2014]のSuppの式(S23)ではExactな値が求まらない(余計に大きい値が求まることがある)
		if (Dm[i][j].first)return Dm[i][j].second;

		int m = -10000000;
		for (int k = i; k <= j; ++k) {
			m = std::max(m, GetDm(i, k - 1) + GetD1(k, j) - C[i][k - 1] - C[k][j]);
			m = std::max(m, GetD1(k, j) - C[k][j]);
		}
		return (Dm[i][j] = std::make_pair(true, m + C[i][j])).second;
	};

	const int ans = GetD(1, n);
	return GetD(1, n);
}

std::vector<std::string>EnumerateStructures(const std::string& sequence, const int max_span, const int max_loop) {
	//可能なRNA二次構造を総当りで全列挙する。

	for (const char c : sequence)assert(c == 'A' || c == 'U' || c == 'G' || c == 'C');

	const int n = int(sequence.size());

	std::function<std::vector<std::string>(std::string)> EnumSubStr;

	EnumSubStr = [&](const std::string subseq) {

		const int nn = int(subseq.size());

		std::vector<std::string>ans;

		if (nn <= TURN + 1) {
			std::string p = "";
			for (int i = 0; i < nn; ++i)p += ".";
			ans.push_back(p);
			return ans;
		}

		//subseq[0]がどれと塩基対を組むかで場合分けする。
		for (int j = TURN + 1; j < nn && j <= max_span; ++j) {
			const int type = parasor_param::GetPairType(subseq[0], subseq[j]);
			if (type == 0)continue;
			const auto x = EnumSubStr(subseq.substr(1, j - 1));
			const auto y = EnumSubStr(subseq.substr(j + 1, nn - j - 1));
			for (const std::string s1 : x)for (const std::string s2 : y) {
				ans.push_back("(" + s1 + ")" + s2);
			}
		}

		//subseq[0]が塩基対を組まない場合
		const auto z = EnumSubStr(subseq.substr(1, nn - 1));
		for (const std::string s : z) {
			ans.push_back("." + s);
		}

		return ans;
	};

	const auto tmp = EnumSubStr(sequence);

	std::vector<std::string>ans;
	for (const std::string s : tmp) {
		if (ComputeMaxLoop(s) > max_loop) {
			continue;
		}
		ans.push_back(s);

		VerifyAndParseStructure(s, sequence, max_span, max_loop);
	}

	return ans;
}

std::vector<std::vector<int>>RandomStructurePK(const std::string& sequence, const int seed) {
	//可能なRNA二次構造のうちどれか一つを返す。
	//返り値は2次元配列で、(i,j)が塩基対を組んでいるなら[i][j]=1で、さもなくば[i][j]=0。i>=jなら常に[i][j]=0。
	//Pseudo-Knotありとする。max_span制約, max_loop制約は考えない。
	//可能な任意の二次構造について、それが返される確率は0より大きいが、必ずしも一様ではない。

	for (const char c : sequence)assert(c == 'A' || c == 'U' || c == 'G' || c == 'C');
	const std::string bp("AU CG GC GU UA UG");

	const int n = int(sequence.size());

	std::vector<std::vector<int>>ans(n + 1, std::vector<int>(n + 1, 0));
	std::vector<int>paired(n, 0);
	std::mt19937_64 rnd(seed);

	for (int i = 0; i < n; ++i)if (paired[i] == 0) {
		std::string p("XX");
		p[0] = sequence[i];
		std::vector<int>candidate;
		candidate.push_back(-1);
		for (int j = i + TURN + 1; j < n; ++j) {
			p[1] = sequence[j];
			if (bp.find(p) != std::string::npos)candidate.push_back(j);
		}
		const int pos = candidate[std::uniform_int_distribution<int>(0, int(candidate.size()) - 1)(rnd)];
		if (pos != -1) {
			ans[i + 1][pos + 1] = 1;
			paired[pos] = 1;
		}
	}

	return ans;
}

int ComputeHammingDistance(const std::string& structure1, const std::string& structure2) {
	//2つのRNA二次構造のハミング距離を返す。

	const int n = int(structure1.size());
	assert(n == int(structure2.size()));

	const std::vector<std::string>str{ structure1,structure2 };
	std::vector<std::set<std::pair<int, int>>>bp(2);
	std::vector<std::stack<int>>pos(2);
	for (int x = 0; x < 2; ++x)for (int i = 0; i < n; ++i) {
		switch (str[x][i]) {
		case '(':
			pos[x].push(i);
			break;
		case ')':
			bp[x].insert(std::make_pair(i, pos[x].top()));
			pos[x].pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}

	std::set<std::pair<int, int>>intersection;
	std::set_intersection(
		bp[0].begin(), bp[0].end(),
		bp[1].begin(), bp[1].end(),
		inserter(intersection, intersection.end()));

	return int(bp[0].size()) + int(bp[1].size()) - 2 * int(intersection.size());
}

int ComputeHammingDistance(const std::vector<std::vector<int>>& structure1, const std::vector<std::vector<int>>& structure2) {
	//2つのRNA二次構造のハミング距離を返す。

	const int n = int(structure1.size());
	assert(n == int(structure2.size()));
	for (int i = 0; i < n; ++i) {
		assert(structure1[i].size() == n);
		assert(structure2[i].size() == n);
	}

	int distance = 0;
	for (int i = 0; i < n; ++i)for (int j = i + 1; j < n; ++j) {
		assert(structure1[i][j] == 0 || structure1[i][j] == 1);
		assert(structure2[i][j] == 0 || structure2[i][j] == 1);
		if (structure1[i][j] != structure2[i][j])++distance;
	}

	return distance;
}

std::string ComputeStructuralContext(const std::string& structure, const int pos) {
	//RNA二次構造structureにおいて、pos番目の塩基(0-origin)のstructural contextを求めて返す。

	const int n = int(structure.size());

	if (structure[pos] != '.')return std::string("Stem");
	if (pos == 0 || pos == n - 1)return std::string("Exterior");

	std::vector<int> bp(n, -1);
	std::stack<int> stk;
	for (int i = 0; i < n; ++i) {
		switch (structure[i]) {
		case '(':
			stk.push(i);
			break;
		case ')':
			bp[i] = stk.top();
			bp[stk.top()] = i;
			stk.pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}

	int dangleflag = 1;
	int danglenum = 0;
	int branchnum = 0;
	for (int i = pos + 1; i != pos; ++i) {
		if (bp[i] != -1) {
			++branchnum;
			danglenum += dangleflag;
			dangleflag = 0;
			i = bp[i];
		}
		else dangleflag = 1;
		if (i == n - 1)return std::string("Exterior");
	}
	if (branchnum == 1)return std::string("Hairpin");
	if (branchnum >= 3)return std::string("Multibranch");
	if (danglenum == 1)return std::string("Bulge");
	return std::string("Interior");
}

double EvalSpecificStructure(const std::string& sequence, const std::string& structure) {
	//特定の配列sequenceが特定の二次構造structureを取るときのボルツマン因子を求めて返す。

	for (const char c : sequence) {
		assert(c == 'A' || c == 'U' || c == 'G' || c == 'C');
	}
	for (const char c : structure) {
		if (!(c == '(' || c == '.' || c == ')')) {
			assert(0);
		}
		assert(c == '(' || c == '.' || c == ')');
	}

	const int n = int(structure.size());

	//塩基対を組んでいる各塩基について、相手を求めてbpに記録する。
	std::vector<int> bp(n, -1);
	std::stack<int> stk;
	for (int i = 0; i < n; ++i) {
		switch (structure[i]) {
		case '(':
			stk.push(i);
			break;
		case ')':
			bp[i] = stk.top();
			bp[stk.top()] = i;
			stk.pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}

	std::function<double(int, int, bool)>factor_in_ij;
	factor_in_ij = [&sequence, &bp, &factor_in_ij](const int i, const int j, const bool ext) {
		//塩基対(i,j)で閉じられる領域のボルツマン因子
		//extは(i,j)が最外側塩基対ならtrueでさもなくばfalse

		//(i,j)で閉じられるループの分岐数を調べて、各分岐の位置を記録する。
		std::vector<std::pair<int, int>>branch;
		for (int p = i + 1; p < j; ++p) {
			if (bp[p] != -1) {
				assert(p < bp[p]);
				branch.push_back(std::make_pair(p, bp[p]));
				p = bp[p];
			}
		}

		double ans = 1.0;
		if (int(branch.size()) == 0) {
			ans *= exp(parasor_param::ParHairpinEnergy(i, j, sequence));
		}
		else if (int(branch.size()) == 1) {
			const int k = branch[0].first;
			const int l = branch[0].second;
			ans *= exp(parasor_param::ParLoopEnergy(i, j, k, l, sequence));
			ans *= factor_in_ij(k, l, false);
		}
		else {
			ans *= exp(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i], sequence[j]), j - 1, i + 1, false, sequence));
			ans *= exp(parasor_param::ParMultiloopClosing());
			int u = j - i - 1;
			for (auto b : branch) {
				const int k = b.first;
				const int l = b.second;
				u -= l - k + 1;
				ans *= exp(parasor_param::ParDangling(parasor_param::GetPairType(sequence[k], sequence[l]), k - 1, l + 1, false, sequence));
				ans *= exp(parasor_param::ParMultiloopInternal());
				ans *= factor_in_ij(k, l, false);
			}
		}

		return ans;
	};

	double ans = 1.0;
	for (int i = 0;; i++) {
		if (bp[i] != -1) {
			assert(i < bp[i]);
			ans *= factor_in_ij(i, bp[i], true);
			ans *= exp(parasor_param::ParDangling(parasor_param::GetPairType(sequence[i], sequence[bp[i]]), i - 1, bp[i] + 1, true, sequence));
			i = bp[i];
		}
		if (i == n - 1)break;
	}

	return ans;
}

double EvalSpecificStructure(const std::string& sequence, const std::vector<std::vector<int>>& structure) {
	//特定の配列sequenceが特定の二次構造structureを取るときのボルツマン因子を求めて返す。
	//引数structureをドット表記に変換して上の関数に入れるよりちょっと速い（はず）。あとMaxLoop制約に違反する場合の処理が違う。

	for (const char c : sequence) {
		assert(c == 'A' || c == 'U' || c == 'G' || c == 'C');
	}

	const int n = int(structure.size());

	assert(structure.size() == n + 1);
	for (int i = 0; i <= n; ++i)assert(structure[i].size() == n + 1);

	//塩基対を組んでいる各塩基について、相手を求めてbpに記録する。
	std::vector<int> bp(n, -1);
	//std::stack<int> stk;
	//for (int i = 0; i < n; ++i) {
	//	switch (structure[i]) {
	//	case '(':
	//		stk.push(i);
	//		break;
	//	case ')':
	//		bp[i] = stk.top();
	//		bp[stk.top()] = i;
	//		stk.pop();
	//		break;
	//	case '.':
	//		break;
	//	default:
	//		assert(0);
	//		break;
	//	}
	//}
	for (int i = 1; i <= n; ++i) {
		int x = -1;
		for (int j = i + 1; j <= n; ++j)if (structure[i][j]) {
			x = j-1;
			break;
		}
		if (x != -1) {
			bp[i - 1] = x;
			bp[x] = i - 1;
		}
	}

	std::function<double(int, int, bool)>factor_in_ij;
	factor_in_ij = [&sequence, &bp, &factor_in_ij](const int i, const int j, const bool ext) {
		//塩基対(i,j)で閉じられる領域のボルツマン因子
		//extは(i,j)が最外側塩基対ならtrueでさもなくばfalse

		//(i,j)で閉じられるループの分岐数を調べて、各分岐の位置を記録する。
		std::vector<std::pair<int, int>>branch;
		for (int p = i + 1; p < j; ++p) {
			if (bp[p] != -1) {
				assert(p < bp[p]);
				branch.push_back(std::make_pair(p, bp[p]));
				p = bp[p];
			}
		}

		double ans = 1.0;
		if (int(branch.size()) == 0) {
			ans *= exp(parasor_param::ParHairpinEnergy(i, j, sequence));
		}
		else if (int(branch.size()) == 1) {
			const int k = branch[0].first;
			const int l = branch[0].second;
			if (k - i + j - l >= 32)return 0.0;//MaxLoop制約に違反する構造入力に対してはassertで落とすのではなくボルツマン因子＝ゼロを返す。
			ans *= exp(parasor_param::ParLoopEnergy(i, j, k, l, sequence));
			ans *= factor_in_ij(k, l, false);
		}
		else {
			ans *= exp(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i], sequence[j]), j - 1, i + 1, false, sequence));
			ans *= exp(parasor_param::ParMultiloopClosing());
			int u = j - i - 1;
			for (auto b : branch) {
				const int k = b.first;
				const int l = b.second;
				u -= l - k + 1;
				ans *= exp(parasor_param::ParDangling(parasor_param::GetPairType(sequence[k], sequence[l]), k - 1, l + 1, false, sequence));
				ans *= exp(parasor_param::ParMultiloopInternal());
				ans *= factor_in_ij(k, l, false);
			}
		}

		return ans;
	};

	double ans = 1.0;
	for (int i = 0;; i++) {
		if (bp[i] != -1) {
			assert(i < bp[i]);
			ans *= factor_in_ij(i, bp[i], true);
			ans *= exp(parasor_param::ParDangling(parasor_param::GetPairType(sequence[i], sequence[bp[i]]), i - 1, bp[i] + 1, true, sequence));
			i = bp[i];
		}
		if (i == n - 1)break;
	}

	return ans;
}

std::string MatrixToDotNotation(const std::vector<std::vector<int>>& structure) {
	//VerifyAndParseStructureの返り値の形式である、バイナリ上三角行列表現のRNA二次構造を引数に取り、
	//ドット表記のRNA二次構造に変換して返す。

	const int N = structure.size() - 1;
	for (int i = 0; i <= N; ++i)assert(structure[i].size() == N + 1);
	for (int i = 0; i <= N; ++i)for (int j = 0; j <= N; ++j) {
		assert(structure[i][j] == 0 || structure[i][j] == 1);
		if (structure[i][j] == 1)assert(1 <= i && i < j);
	}

	std::string answer;
	for (int i = 0; i < N; ++i)answer += std::string(".");
	for (int i = 0; i < N; ++i)for (int j = 0; j < N; ++j)if (structure[i + 1][j + 1] == 1) {
		assert(answer[i] == '.' && answer[j] == '.');
		answer[i] = '(';
		answer[j] = ')';
	}

	return answer;
}

}
