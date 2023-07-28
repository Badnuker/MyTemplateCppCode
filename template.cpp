#include<bits/stdc++.h>
#define int long long
using namespace std;
const int maxn = 1e5 + 5;

namespace bignum {
	class Big {
		private:
			int len, num[maxn];

		public:
			Big(int x = 0) {
				len = 0;
				memset(num, 0, sizeof(num));
				while (x) {
					num[++len] = x % 10;
					x /= 10;
				}
			}

			void read() {
				string str;
				cin >> str;
				len = str.length();
				reverse(str.begin(), str.end());
				str = ' ' + str;
				for (int i = 1; i <= len; i++) {
					num[i] = str[i] - '0';
				}
			}

			void print() {
				for (int i = len; i >= 1; i--) {
					cout << num[i];
				}
			}

			Big operator+(const Big &x)const {
				Big res;
				res.len = max(len, x.len);
				for (int i = 1; i <= res.len; i++) {
					res.num[i] = num[i] + x.num[i];
				}
				for (int i = 1; i <= res.len; i++) {
					res.num[i + 1] += res.num[i] / 10;
					res.num[i] %= 10;
				}
				while (res.num[res.len] == 0 && res.len > 1) {
					res.len--;
				}
				return res;
			}

			Big operator -(const Big &x)const {
				Big tmp;
				tmp.len = max(len, x.len);
				for (int i = 1; i <= tmp.len; i++) {
					tmp.num[i] = num[i] - x.num[i];
				}
				for (int i = 1; i < tmp.len; i++) {
					if (tmp.num[i] < 0) {
						tmp.num[i] += 10;
						tmp.num[i + 1]--;
					}
				}
				while (tmp.num[tmp.len] == 0) {
					tmp.len--;
				}
				return tmp;
			}

			Big operator*(const Big &x)const {
				Big res;
				res.len = len + x.len;
				for (int i = 1; i <= len; i++) {
					for (int j = 1; j <= x.len; j++) {
						res.num[i + j - 1] += num[i] * x.num[j];
					}
				}
				for (int i = 1; i <= res.len; i++) {
					res.num[i + 1] += res.num[i] / 10;
					res.num[i] %= 10;
				}
				while (res.num[res.len] == 0 && res.len > 1) {
					res.len--;
				}
				return res;
			}

			Big operator+(const int &x)const {
				return (*this) + Big(x);
			}

			Big operator-(const int &x)const {
				return (*this) - Big(x);
			}

			Big operator*(const int &x)const {
				return (*this) * Big(x);
			}

			bool operator==(const Big &x)const {
				if (len != x.len) {
					return 0;
				}
				for (int i = len; i; i--) {
					if (num[i] != x.num[i]) {
						return 0;
					}
				}
				return 1;
			}

			bool operator<(const Big &x)const {
				if (len < x.len) {
					return 1;
				}
				if (len > x.len) {
					return 0;
				}
				for (int i = len; i; i--) {
					if (num[i] < x.num[i]) {
						return 1;
					}
					if (num[i] > x.num[i]) {
						return 0;
					}
				}
				return 0;
			}

			bool operator >(const Big &x)const {
				return !((*this) < x || (*this) == x);
			}

			bool operator <=(const Big &x)const {
				return !((*this) > x);
			}

			bool operator >=(const Big &x)const {
				return !((*this) < x);
			}

	};

}

namespace math {
	int lowbit(int x) {
		return x & (-x);
	}

	int ksm(int x, int y) {
		int res = 1;
		while (y) {
			if (y & 1) {
				res *= x;
			}
			x *= x;
			y >>= 1;
		}
		return res;
	}

}

namespace bittree {
	int n, m;

	class BitTree {
		private:
			int c[maxn];

		public:
			void update(int pos, int val) {
				for (; pos <= n; pos += math::lowbit(pos)) {
					c[pos] += val;
				}
			}

			int query(int pos) {
				int res = 0;
				for (; pos; pos -= math::lowbit(pos)) {
					res += c[pos];
				}
				return res;
			}

	} bit;

}

namespace segtree {
#define lp (p<<1)
#define rp (p<<1|1)
#define mid ((l+r)>>1)

	int n, m;
	int a[maxn];

	class AddSumSegTree {
		private:
			int sum[maxn << 2], tag[maxn << 2];

			void pushup(int p) {
				sum[p] = sum[lp] + sum[rp];
			}

			void pushdown(int l, int r, int p) {
				if (tag[p]) {
					tag[lp] += tag[p];
					tag[rp] += tag[p];
					sum[lp] += (mid - l + 1) * tag[p];
					sum[rp] += (r - mid) * tag[p];
					tag[p] = 0;
				}
			}
		public:
			void build(int l = 1, int r = n, int p = 1) {
				if (l == r) {
					sum[p] = a[l];
					return;
				}
				build(l, mid, lp);
				build(mid + 1, r, rp);
				pushup(p);
			}

			void update(int L, int R, int c, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					sum[p] += (r - l + 1) * c;
					tag[p] += c;
					return;
				}
				pushdown(l, r, p);
				if (L <= mid) {
					update(L, R, c, l, mid, lp);
				}
				if (R > mid) {
					update(L, R, c, mid + 1, r, rp);
				}
				pushup(p);
			}

			int query(int L, int R, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					return sum[p];
				}
				pushdown(l, r, p);
				int res = 0;
				if (L <= mid) {
					res += query(L, R, l, mid, lp);
				}
				if (R > mid) {
					res += query(L, R, mid + 1, r, rp);
				}
				return res;
			}
	};

	class ModSumSegTree {
		private:
			int sum[maxn << 2], tag[maxn << 2];

			void pushup(int p) {
				sum[p] = sum[lp] + sum[rp];
			}

			void pushdown(int l, int r, int p) {
				if (tag[p] != -1) {
					tag[lp] = tag[p];
					tag[rp] = tag[p];
					sum[lp] = (mid - l + 1) * tag[p];
					sum[rp] = (r - mid) * tag[p];
					tag[p] = -1;
				}
			}
		public:
			void build(int l = 1, int r = n, int p = 1) {
				if (l == r) {
					sum[p] = a[l];
					return;
				}
				build(l, mid, lp);
				build(mid + 1, r, rp);
				pushup(p);
			}

			void update(int L, int R, int c, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					sum[p] = (r - l + 1) * c;
					tag[p] = c;
					return;
				}
				pushdown(l, r, p);
				if (L <= mid) {
					update(L, R, c, l, mid, lp);
				}
				if (R > mid) {
					update(L, R, c, mid + 1, r, rp);
				}
				pushup(p);
			}

			int query(int L, int R, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					return sum[p];
				}
				pushdown(l, r, p);
				int res = 0;
				if (L <= mid) {
					res += query(L, R, l, mid, lp);
				}
				if (R > mid) {
					res += query(L, R, mid + 1, r, rp);
				}
				return res;
			}
	};
#undef lp
#undef rp
#undef mid

}

namespace heavychain {
#define lp (p<<1)
#define rp (p<<1|1)
#define mid ((l+r)>>1)

	int n, m, dfn_id;
	int a[maxn];
	int fa[maxn], sz[maxn], dep[maxn];
	int son[maxn], top[maxn], dfn[maxn], rnk[maxn];
	vector<int> p[maxn];

	void dfs1(int u) {
		sz[u] = 1;
		for (auto v : p[u]) {
			if (v == fa[u]) {
				continue;
			}
			fa[v] = u;
			dep[v] = dep[u] + 1;
			dfs1(v);
			sz[u] += sz[v];
			if (sz[v] > sz[son[u]]) {
				son[u] = v;
			}
		}
	}

	void dfs2(int u, int t) {
		top[u] = t;
		dfn[u] = ++dfn_id;
		rnk[dfn_id] = u;
		if (!son[u]) {
			return;
		}
		dfs2(son[u], t);
		for (auto v : p[u]) {
			if (v == fa[u] || v == son[u]) {
				continue;
			}
			dfs2(v, v);
		}
	}


	class AddSumSegTree {
		private:
			int sum[maxn << 2], tag[maxn << 2];

			void pushup(int p) {
				sum[p] = sum[lp] + sum[rp];
			}

			void pushdown(int l, int r, int p) {
				if (tag[p]) {
					tag[lp] += tag[p];
					tag[rp] += tag[p];
					sum[lp] += (mid - l + 1) * tag[p];
					sum[rp] += (r - mid) * tag[p];
					tag[p] = 0;
				}
			}

		public:
			void build(int l = 1, int r = n, int p = 1) {
				if (l == r) {
					sum[p] = a[rnk[l]];
					return;
				}
				build(l, mid, lp);
				build(mid + 1, r, rp);
				pushup(p);
			}

			void update(int L, int R, int c, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					sum[p] += (r - l + 1) * c;
					tag[p] += c;
					return;
				}
				pushdown(l, r, p);
				if (L <= mid) {
					update(L, R, c, l, mid, lp);
				}
				if (R > mid) {
					update(L, R, c, mid + 1, r, rp);
				}
				pushup(p);
			}

			int query(int L, int R, int l = 1, int r = n, int p = 1) {
				if (L <= l && r <= R) {
					return sum[p];
				}
				pushdown(l, r, p);
				int res = 0;
				if (L <= mid) {
					res += query(L, R, l, mid, lp);
				}
				if (R > mid) {
					res += query(L, R, mid + 1, r, rp);
				}
				return res;
			}
	} asst;
#undef lp
#undef rp
#undef mid

	void AddChain(int x, int y, int z) {
		while (top[x] != top[y]) {
			if (dep[top[x]] < dep[top[y]]) {
				swap(x, y);
			}
			asst.update(dfn[top[x]], dfn[x], z, 1);
			x = fa[top[x]];
		}
		if (dep[x] > dep[y]) {
			swap(x, y);
		}
		asst.update(dfn[x], dfn[y], z, 1);
	}

	int QueryChainSum(int x, int y) {
		int res = 0;
		while (top[x] != top[y]) {
			if (dep[top[x]] < dep[top[y]]) {
				swap(x, y);
			}
			res += asst.query(dfn[top[x]], dfn[x], 1);
			x = fa[top[x]];
		}
		if (dep[x] > dep[y]) {
			swap(x, y);
		}
		res += asst.query(dfn[x], dfn[y], 1);
		return res;
	}

}

using namespace heavychain;

signed main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	return 0;
}
