! function e(t, n, r) {
  function i(a, s) {
    if (!n[a]) {
      if (!t[a]) {
        var u = "function" == typeof require && require;
        if (!s && u) return u(a, !0);
        if (o) return o(a, !0);
        var l = new Error("Cannot find module '" + a + "'");
        throw l.code = "MODULE_NOT_FOUND", l
      }
      var c = n[a] = {
        exports: {}
      };
      t[a][0].call(c.exports, (function(e) {
        return i(t[a][1][e] || e)
      }), c, c.exports, e, t, n, r)
    }
    return n[a].exports
  }
  for (var o = "function" == typeof require && require, a = 0; a < r.length; a++) i(r[a]);
  return i
}({
  1: [function(e, t, n) {
    "use strict";
    t.exports = function(e, t) {
      e.html5Mode(!1), e.hashPrefix("!"), t.debugInfoEnabled(!1)
    }
  }, {}],
  2: [function(e, t, n) {
    "use strict";
    e("jquery");
    var r = e("angular");
    e("angular-route"), e("angular-sanitize"), e("angular-animate"), e("angular-ui-bootstrap");
    window._ = e("underscore"), r.module("underscore", []).factory("_", ["$window", function(e) {
      return e._
    }]);
    var i = r.module("tmcApp", ["ngRoute", "ngSanitize", "ui.bootstrap", "underscore", "ngAnimate"]);
    e("./services"), e("./controllers"), e("./directives"), i.config(["$routeProvider", e("./app.route")]), i.config(["$locationProvider",
      "$compileProvider", e("./app.config")
    ]), i.run(["$rootScope", "$uibModal", "modalService", e("./app.run")])
  }, {
    "./app.config": 1,
    "./app.route": 3,
    "./app.run": 4,
    "./controllers": 7,
    "./directives": 11,
    "./services": 18,
    angular: 32,
    "angular-animate": 24,
    "angular-route": 26,
    "angular-sanitize": 28,
    "angular-ui-bootstrap": 30,
    jquery: 33,
    underscore: 34
  }],
  3: [function(e, t, n) {
    "use strict";
    t.exports = function(e) {
      e.when("/", {
        redirectTo: "/main"
      }).when("/main", {
        templateUrl: "views/main-b768bd1c14.html",
        controller: "MainCtrl",
        resolve: {
          prodData: ["tmData", function(e) {
            return e.getData("data/tmcalculatordata.json")
          }]
        }
      }).when("/batch", {
        templateUrl: "views/batch-fe85a81f13.html",
        controller: "BatchCtrl",
        resolve: {
          prodData: ["tmData", function(e) {
            return e.getData("data/tmcalculatordata.json")
          }]
        }
      }).when("/testbs", {
        templateUrl: "views/testbs.html",
        controller: "TestbsCtrl"
      }).when("/about", {
        templateUrl: "views/about-f46b5427a5.html",
        controller: "AboutCtrl"
      }).when("/help", {
        templateUrl: "views/help-8ae8e24615.html",
        controller: "AboutCtrl"
      }).when("/history", {
        templateUrl: "views/history-56a6a55424.html",
        controller: "AboutCtrl"
      }).when("/:p1seq/:p2seq", {
        templateUrl: "views/main-b768bd1c14.html",
        controller: "MainCtrl",
        resolve: {
          prodData: ["tmData", function(e) {
            return e.getData("data/tmcalculatordata.json")
          }]
        }
      }).when("/:p1seq", {
        templateUrl: "views/main-b768bd1c14.html",
        controller: "MainCtrl",
        resolve: {
          prodData: ["tmData", function(e) {
            return e.getData("data/tmcalculatordata.json")
          }]
        }
      }).otherwise({
        redirectTo: "/"
      })
    }
  }, {}],
  4: [function(e, t, n) {
    "use strict";
    var r = e("jquery");
    t.exports = function(e, t, n) {
      r.event.addProp("dataTransfer"), e.appVersion = "1.16.5", e.state = {}, e.state.input = {}, e.state.result = {}, e.state.results1 = [], e.state
        .results2 = [], e.state.input.p1 = "", e.state.input.p2 = "", e.state.input.ct = .25, e.state.input.groupKeys = [], e.state.input.products = {},
        e.state.input.prodKeys = [], e.state.input.group = "", e.state.input.product = "", e.state.result.tm1 = "---", e.state.result.tm2 = "---", e
        .state.result.len1 = "---", e.state.result.len2 = "---", e.state.result.gc1 = "---", e.state.result.gc2 = "---", e.state.result.ta = "---", e
        .state.result.itemlist = [], e.state.result.critlist = [], e.state.p1status = "", e.state.p2status = "", e.state.results1.showexpanded1 = !1, e
        .state.results2.showexpanded1 = !1;
      var i = window,
        o = document,
        a = o.documentElement,
        s = o.getElementsByTagName("body")[0],
        u = i.innerWidth || a.clientWidth || s.clientWidth;
      i.innerHeight || a.clientHeight || s.clientHeight;
      e.initScreenWidth = u, e.screenMessage = !1, e.openmodal = function(e, n) {
        var r = {
          templateUrl: "views/modals/" + e + ".html",
          size: n
        };
        return t.open(r, {})
      }, e.$watch((function() {
        return MathJax.Hub.Queue(["Typeset", MathJax.Hub]), !0
      }))
    }
  }, {
    jquery: 33
  }],
  5: [function(e, t, n) {
    "use strict";
    e("angular");
    t.exports = function(e) {}
  }, {
    angular: 32
  }],
  6: [function(e, t, n) {
    "use strict";
    e("angular");
    var r = e("../vendorlibs/FileSaver");
    t.exports = function(e, t, n, i, o, a, s, u, l) {
      e.input = {}, e.result = {}, e.output = [], e.input.p1 = "", e.input.p2 = "", e.input.id1 = "", e.input.id2 = "", e.input.ct = .25, e.result.tm1 =
        "", e.result.tm2 = "", e.result.ta = "", e.result.itemlist = [], e.result.critlist = [], e.input.groupKeys = [], e.input.products = {}, e.input
        .prodKeys = [], e.input.group = "", e.input.product = "", e.input.batch = "", e.input.interleaved = !1, e.input.filename = "", e.result.batch =
        "", e.result.batch2 = "", e.p1status = "", e.p2status = "", console.log(n), i.setData(n), console.log(i.tmcdata), i.restoreUserPrefs() ? (e
          .lastp1 = i.restoreUserPrefs().p1, e.lastp2 = i.restoreUserPrefs().p2) : (e.lastp1 = "", e.lastp2 = "");
      var c = s.p1seq || e.lastp1 || "",
        d = s.p2seq || e.lastp2 || "";
      e.prefill = function() {
        e.input.interleaved ? e.input.batch =
          "P1fwd; AGCGGATAACAATTTCACACAGGA\nP1rev; GTAAAACGACGGCCAGT\nP3fwd; AGCGGATAAGGGCAATTTCAC\nP3rev; GTAAAACGACGGCCA\n" : e.input.batch =
          "P1fwd; AGCGGATAACAATTTCACACAGGA; P1rev; GTAAAACGACGGCCAGT\npBeta; AGCGGATAACAATTTCAC\nP3fwd; AGCGGATAAGGGCAATTTCAC; P3rev; GTAAAACGACGGCCA\n",
          e.output.showresultstable = !1, e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show", e.runCalc3(e.input.interleaved)
      }, e.clearCalc = function() {
        e.input.batch = "", e.input.filename = "", e.output = [], e.result.itemlist = [], e.result.critlist = [], e.result.batch = "", e.result
          .batch2 = "", e.runmsg = "", document.getElementById("fileinput").value = "", e.output.showresultstable = !1, e.data_toggle_label = e.output
          .showresultstable ? "Hide" : "Show"
      }, e.setGroups = function() {
        e.input.groupKeys = i.getGroupKeys(), e.input.group = e.input.groupKeys[0]
      }, e.setProducts = function() {
        e.input.group;
        e.input.productKeys = i.getProductKeysForGroup(e.input.group), e.input.products = i.getProductsForGroup(e.input.group), e.input.product = e
          .input.products[0].id;
        e.input.product;
        e.setCt()
      }, e.setCt = function() {
        var t = e.input.group,
          n = e.input.product;
        0 === t.indexOf("Phusion") || 0 === t.indexOf("Q5") ? e.input.ct = 500 : 0 === t.indexOf("LongAmp") ? e.input.ct = 400 : 0 === t.indexOf(
          "Master") ? 0 === n.indexOf("phusion") || 0 === n.indexOf("q5") ? e.input.ct = 500 : 0 === n.indexOf("lataq") || 0 === n.indexOf(
          "lahstaq") ? e.input.ct = 400 : e.input.ct = 200 : e.input.ct = 200, e.runCalc3(e.input.interleaved)
      }, e.runCalc3 = function(t) {
        t = !!t, console.log("Interleaved", t);
        for (var n, r, s, u, l, c, d, p, f, h, g, m, v, $, b, y, w, x, C, S = "", A = "", T = "", k = "", D = "", M = e.input.ct / 1e3, E = e.input
            .product, O = e.input.group, N = [], P = [], q = 0, I = null, L = 1, R = e.input.batch, j = R.split(/\n\r?/), _ = 0; _ < j.length; ++_) j[
          _].match(/^\s*$/) ? j.splice(_, 1) : j[_] = j[_].replace(/[;,\t\s]$/, "");
        if (j && 0 !== j.length) {
          console.log("after cleanup");
          for (var V = j[0].slice(), U = V.split(""), H = 0; H < U.length; ++H) console.log(H, U[H], V.charCodeAt(H));
          var F = {};
          if (console.log("testing separators"), F.semicolon = V.split(";").length, console.log("fields when split with semicolon: ", F.semicolon), F
            .comma = V.split(",").length, console.log("fields when split with comma: ", F.comma), F.tab = V.split("\t").length, console.log(
              "fields when split with tab: ", F.tab), F.ws = V.split(/[\s\W]+/).length, console.log("fields when split with ws: ", F.ws), t ? (I =
              2 === F.semicolon ? ";" : 2 === F.comma ? "," : 2 === F.tab ? "\t" : 2 === F.ws ? /[\s\W]+/ : null, L = 2) : (I = 2 === F.semicolon ||
              3 === F.semicolon || 4 === F.semicolon ? ";" : 2 === F.comma || 3 === F.comma || 4 === F.comma ? "," : 2 === F.tab || 3 === F.tab ||
              4 === F.tab ? "\t" : 2 === F.ws || 3 === F.ws || 4 === F.ws ? /[\s\W]+/ : null, L = 1), null == I) return P = [
              "Unable to read primers and/or invalid data format. "
            ], e.result.critlist = P, e.result.itemlist = [], R = "", e.input.filename = "", document.getElementById("fileinput").value = "", e
            .runmsg = "", void(e.output = []);
          if (e.result = {}, e.output = [], e.ctstatus = "", e.runmsg = "", e.output.showresultstable = !1, e.data_toggle_label = e.output
            .showresultstable ? "Hide" : "Show", M <= 0 || isNaN(M / 1)) return !1, e.ctstatus = "invalidct", P.push(
            "Invalid primer concentration. "), e.result.itemlist = N, e.result.critlist = P, void(e.p1status = "");
          switch ("onetaq_gc" === (s = i.getBufferIdForProduct(E)) && (q = 5), u = i.getBufferSaltForProduct(e.input.product), O) {
            case "Phusion":
            case "Phusion Hot Start Flex":
              l = 5, M /= 4;
              break;
            case "Master Mix":
              0 === E.indexOf("phusion") ? (l = 5, M /= 4) : l = 4;
              break;
            default:
              l = 4
          }($ = i.validateBuffer(E, O, s)).hasCritWarnings && Array.prototype.push.apply(P, $.critwarnings), $.hasWarnings && Array.prototype.push
            .apply(N, $.warnings), e.result.itemlist = N, e.result.critlist = P, o.setCt(M).setMonosalt(u).setMethod(l).setDMSO(q), d = 0, p = 0, g =
            0;
          for (var B = 0; B <= j.length - L; B += L) {
            if (D = "OK", f = !1, h = !1, t) {
              if (n = j[B].split(I), r = j[B + 1].split(I), n.length < 2 && r.length < 2) {
                ++g;
                continue
              }
              2 === n.length ? (S = n[1], T = n[0]) : (S = "", T = ""), 2 === r.length ? (A = r[1], k = r[0]) : (A = "", k = "")
            } else if (4 === (n = j[B].split(I)).length) S = n[1].replace(/[/W/s]/g, ""), A = n[3].replace(/[/W/s]/g, ""), T = n[0], k = n[2];
            else if (3 === n.length) S = n[1].replace(/[/W/s]/g, ""), A = n[2].replace(/[/W/s]/g, ""), T = n[0], k = n[0];
            else {
              if (2 !== n.length) {
                ++g;
                continue
              }
              S = n[1].replace(/[/W/s]/g, ""), A = "", T = n[0], k = ""
            }
            if (e.input.p1 = S.toUpperCase(), e.input.p2 = A.toUpperCase(), e.input.id1 = T, e.input.id2 = k, (m = i.validateInput(e.input.p1.replace(
                /\s/g, ""), e.input.p2.replace(/\s/g, ""), M, E, O)).p1isValid, m.p2isValid, c = m.hasCritWarnings, e.p1status = m.p1status, e
              .p2status = m.p2status, c && (D = m.critwarnings.join("-- "), f = !0), "invalidseq" !== e.p1status) {
              var G, W, z = a.expandDNASeq(S.replace(/\s/g, ""));
              if (1 === z.length) y = (W = o.setCt(M).setMonosalt(u).setMethod(l).setDMSO(q).calcTm_nomm(S.replace(/\s/g, ""))).tm, y = Math.round(
                Math.round(10 * y) / 10), e.result.tm1 = y, e.result.len1 = W.len, e.result.gc1 = 100 * W.fgc, e.result.gc1 = Math.round(10 * e
                .result.gc1 / 10);
              else {
                e.results1 = [], o.setCt(M).setMonosalt(u).setMethod(l).setDMSO(q);
                for (var Y = 0; Y < z.length; ++Y) G = z[Y], W = o.calcTm_nomm(G), y = Math.round(Math.round(10 * W.tm) / 10), e.results1.push({
                  seq: G,
                  tm: y,
                  len: G.length,
                  gc: Math.round(100 * W.fgc * 10 / 10)
                });
                y = e.results1.reduce((function(e, t) {
                    return t.tm < e && (e = t.tm), e
                  }), e.results1[0].tm), x = e.results1.reduce((function(e, t) {
                    return t.tm > e && (e = t.tm), e
                  }), e.results1[0].tm), e.result.tm1 = y !== x ? y + " - " + x : y, e.result.len1 = S.replace(/\s/g, "").length, e.result.gc1 =
                  "---", e.result.tm1m = x
              }
            }
            if ("invalidseq" !== e.p2status) {
              var K, Q, X = a.expandDNASeq(A.replace(/\s/g, ""));
              if (1 === X.length) w = (Q = o.setCt(M).setMonosalt(u).setMethod(l).setDMSO(q).calcTm_nomm(A.replace(/\s/g, ""))).tm, w = Math.round(
                Math.round(10 * w) / 10), e.result.tm2 = w, e.result.len2 = Q.len, e.result.gc2 = 100 * Q.fgc, e.result.gc2 = Math.round(10 * e
                .result.gc2 / 10);
              else {
                e.results2 = [], o.setCt(M).setMonosalt(u).setMethod(l).setDMSO(q);
                for (var J = 0; J < z.length; ++J) K = X[J], Q = o.calcTm_nomm(K), w = Math.round(Math.round(10 * Q.tm) / 10), e.results2.push({
                  seq: K,
                  tm: w,
                  len: K.length,
                  gc: Math.round(100 * Q.fgc * 10 / 10)
                });
                w = e.results2.reduce((function(e, t) {
                    return t.tm < e && (e = t.tm), e
                  }), e.results2[0].tm), C = e.results2.reduce((function(e, t) {
                    return t.tm > e && (e = t.tm), e
                  }), e.results2[0].tm), e.result.tm2 = w !== C ? w + " - " + C : w, e.result.len2 = A.replace(/\s/g, "").length, e.result.gc2 =
                  "---", e.result.tm2m = C
              }
            }
            "invalidseq" !== e.p1status && "invalidseq" !== e.p2status ? (b = i.getAnnealTemp(S, y, A, w, O, E), e.result.ta = Math.round(b), (v = i
                .validateTm(y, w, b, E, O, s)).hasCritWarnings && (h = !0, D += "--" + v.critwarnings.join("--")), v.hasWarnings && (h = !0, D +=
                "-- " + v.warnings.join("-- "))) : e.result.ta = "---", f && (d += 1), h && (p += 1), "invalidseq" !== e.p1status && e.output.push([e
                .input.id1, e.input.p1, e.result.tm1, e.result.ta, D
              ].join("\t")), "invalidseq" !== e.p2status && e.output.push([e.input.id2, e.input.p2, e.result.tm2, e.result.ta, D].join("\t")), e
              .p1status = ""
          }
          e.output.length, e.runmsg = g > 0 ? g + " Invalid line(s) " : "", e.runmsg += e.output.length + " primers processed. Errors: " + d +
            " Warnings: " + p, e.result.batch = e.output.join("\n"), e.result.batch2 = [];
          for (var Z = 0; Z < e.output.length; ++Z) e.result.batch2.push(e.output[Z].split("\t"));
          0 == e.output.length && j.length > 0 ? e.novaliddatamsg = "Unable to detect format or no valid data entered -- check format of 1st line." :
            e.novaliddatamsg = ""
        }
      }, e.getBatchResults = function() {
        return e.result.batch
      }, e.downloadData = function() {
        var t = e.getBatchResults(),
          n = new Blob([t], {
            type: "text/plain;charset=utf-8"
          });
        r.saveAs(n, "tmcalc_batch.txt")
      };
      e.$on("$routeChangeStart", (function() {})), e.about = function() {
          u.path("#/about")
        }, e.safeApply = function(e) {
          var t = this.$root.$$phase;
          "$apply" == t || "$digest" == t ? e && "function" == typeof e && e() : this.$apply(e)
        }, e.processFileData = function(t) {
          e.safeApply((function() {
            e.input.batch = t, console.log("filedata:" + t), e.output.showresultstable = !1, e.data_toggle_label = e.output.showresultstable ?
              "Hide" : "Show", e.runCalc3(e.input.interleaved)
          }))
        }, e.switch2single = function() {
          u.path("/")
        }, e.toggle_data_display = function() {
          e.output.showresultstable = !e.output.showresultstable, e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show"
        }, e.output.showresultstable = !1, e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show", e.input.p1 = c, e.input.p2 = d, e
        .setGroups(), e.setProducts()
    }
  }, {
    "../vendorlibs/FileSaver": 22,
    angular: 32
  }],
  7: [function(e, t, n) {
    "use strict";
    var r = e("angular").module("tmcApp");
    r.controller("MainCtrl", ["$scope", "$rootScope", "nebutil", "prodData", "tmcalculatorData", "tmcalc", "neb_bioseq", "$routeParams", "$location",
      "$window", e("./main")
    ]), r.controller("BatchCtrl", ["$scope", "nebutil", "prodData", "tmcalculatorData", "tmcalc", "neb_bioseq", "$routeParams", "$location", "$window",
      e("./batch")
    ]), r.controller("AboutCtrl", ["$scope", e("./about")])
  }, {
    "./about": 5,
    "./batch": 6,
    "./main": 8,
    angular: 32
  }],
  8: [function(e, t, n) {
    "use strict";
    e("angular");
    t.exports = function(e, t, n, r, i, o, a, s, u, l) {
      e.input = {}, e.result = {}, e.results1 = [], e.results2 = [], e.input.p1 = "", e.input.p2 = "", e.input.ct = .25, e.result.tm1 = "---", e.result
        .tm2 = "---", e.result.len1 = "---", e.result.len2 = "---", e.result.gc1 = "---", e.result.gc2 = "---", e.result.ta = "---", e.result
        .itemlist = [], e.result.critlist = [], e.input.groupKeys = [], e.input.products = {}, e.input.prodKeys = [], e.input.group = "", e.input
        .product = "", e.p1status = "", e.p2status = "", console.log(r), i.setData(r), console.log(i.tmcdata), i.restoreUserPrefs() ? (e.lastp1 = i
          .restoreUserPrefs().p1, e.lastp2 = i.restoreUserPrefs().p2) : (e.lastp1 = "", e.lastp2 = "");
      var c = s.p1seq || e.lastp1 || "";
      s.p2seq || e.lastp2;
      console.log("p1seq " + c), e.getUrl = function(e) {
        var t = getPolForName(e),
          n = null;
        return t && (n = t.URL), console.log("URL " + n), n
      };
      e.prefill = function() {
        console.log("prefill"), e.input.p1 = "AGCGGATAACAATTTCACACAGGA", e.input.p2 = "GTA AAA CGA CGG CCA GT", e.runCalc2()
      }, e.clearCalc = function() {
        console.log("clearCalc"), e.input.p1 = "", e.input.p2 = "", e.result.tm1 = "---", e.result.tm2 = "---", e.result.len1 = "---", e.result.len2 =
          "---", e.result.gc1 = "---", e.result.gc2 = "---", e.result.ta = "---", e.result.itemlist = [], e.result.critlist = [], e.results1 = [], e
          .results2 = [], e.p1status = "", e.p2status = "", e.showInfoMessage = !1
      }, e.setGroups = function() {
        console.log("setGroups called"), e.input.groupKeys = i.getGroupKeys(), e.input.group = e.input.groupKeys[0]
      }, e.setProducts = function() {
        var t;
        console.log("setProducts called"), t = e.input.group, e.input.productKeys = i.getProductKeysForGroup(t), e.input.products = i
          .getProductsForGroup(t), e.input.product = e.input.products[0].id, e.input.product, e.setCt()
      }, e.setCt = function() {
        console.log("setCt called");
        var t = e.input.group,
          n = e.input.product;
        0 === t.indexOf("Phusion") || 0 === t.indexOf("Q5") ? e.input.ct = 500 : 0 === t.indexOf("LongAmp") ? e.input.ct = 400 : 0 === t.indexOf(
          "Master") ? 0 === n.indexOf("phusion") || 0 === n.indexOf("q5") ? e.input.ct = 500 : 0 === n.indexOf("lataq") || 0 === n.indexOf(
          "lahstaq") ? e.input.ct = 400 : e.input.ct = 200 : e.input.ct = 200, e.runCalc2()
      }, e.runCalc2 = function() {
        console.log("runCalc2");
        var t, r, s, u, l, c, d, p, f, h, g = e.input.p1,
          m = e.input.p2,
          v = e.input.ct / 1e3,
          $ = e.input.product,
          b = e.input.group,
          y = n.isNumeric(v),
          w = [],
          x = [],
          C = 0,
          S = n.isValidSeq(g, !1),
          A = n.isValidSeq(m, !1),
          T = i.validateInput(g.replace(/\s/g, ""), m.replace(/\s/g, ""), v, $, b);
        if (y = T.ctisValid, S = T.p1isValid, A = T.p2isValid, T.hasCritWarnings, T.hasWarnings, x = T.critwarnings, w = T.warnings, e.p1status = T
          .p1status, e.p2status = T.p2status, e.result.itemlist = w, e.result.critlist = x, e.ctstatus = "", e.results1 = [], e.results2 = [], e
          .results2.showexpanded2 = !1, e.results1.showexpanded1 = !1, e.expanded2_toggle_label = e.results2.showexpanded2 ? "Hide" : "Show", e
          .expanded1_toggle_label = e.results1.showexpanded1 ? "Hide" : "Show", !y) return e.result.ta = "---", e.result.tm1 = null, e.result.tm2 =
          null, e.ctstatus = "invalidct", S && (e.p1status = ""), void(A && (e.p2status = ""));
        if (console.log("checkpoint 1"), console.log("p1 is: " + g), console.log("p2 is: " + m), console.log("p1 is valid: " + S), console.log(
            "p2 is valid: " + A), console.log("p1status: " + e.p1status), console.log("p2status: " + e.p2status), "invalidseq" === e.p1status &&
          "invalidseq" === e.p2status) return e.result.tm1 = "---", e.result.tm2 = "---", e.result.len1 = "---", e.result.len2 = "---", e.result.gc1 =
          "---", e.result.gc2 = "---", e.result.ta = "---", S && 0 === g.length && (e.p1status = ""), A && 0 === m.length && (e.p2status = ""),
          void console.log("checkpoint 1a");
        switch ("invalidseq" === e.p1status && (e.result.tm1 = "---", e.result.len1 = "---", e.result.gc1 = "---", e.result.ta = "---"),
          "invalidseq" === e.p2status && (e.result.tm2 = "---", e.result.len2 = "---", e.result.gc2 = "---", e.result.ta = "---"), console.log(
            "checkpoint 2"), "onetaq_gc" === (t = i.getBufferIdForProduct($)) && (C = 5), r = i.getBufferSaltForProduct(e.input.product), b) {
          case "Phusion":
          case "Phusion Hot Start Flex":
            s = 5, v /= 4;
            break;
          case "Master Mix":
            0 === $.indexOf("phusion") ? (s = 5, v /= 4) : s = 4;
            break;
          default:
            s = 4
        }
        if ("invalidseq" !== e.p1status) {
          var k, D = a.expandDNASeq(g.replace(/\s/g, ""));
          if (1 === D.length) l = (f = o.setCt(v).setMonosalt(r).setMethod(s).setDMSO(C).calcTm_nomm(g.replace(/\s/g, ""))).tm, l = Math.round(Math
            .round(10 * l) / 10), e.result.tm1 = l, e.result.len1 = f.len, e.result.gc1 = 100 * f.fgc, e.result.gc1 = Math.round(10 * e.result.gc1 /
            10);
          else {
            e.results1 = [], o.setCt(v).setMonosalt(r).setMethod(s).setDMSO(C);
            for (var M = 0; M < D.length; ++M) k = D[M], f = o.calcTm_nomm(k), l = Math.round(Math.round(10 * f.tm) / 10), e.results1.push({
              seq: k,
              tm: l,
              len: k.length,
              gc: Math.round(100 * f.fgc * 10 / 10)
            });
            l = e.results1.reduce((function(e, t) {
                return t.tm < e && (e = t.tm), e
              }), e.results1[0].tm), d = e.results1.reduce((function(e, t) {
                return t.tm > e && (e = t.tm), e
              }), e.results1[0].tm), e.result.tm1 = l !== d ? l + " - " + d : l, e.result.len1 = g.replace(/\s/g, "").length, e.result.gc1 = "---", e
              .result.tm1m = d
          }
        }
        if ("invalidseq" !== e.p2status) {
          var E, O = a.expandDNASeq(m.replace(/\s/g, ""));
          if (1 === O.length) c = (h = o.setCt(v).setMonosalt(r).setMethod(s).setDMSO(C).calcTm_nomm(m.replace(/\s/g, ""))).tm, c = Math.round(Math
            .round(10 * c) / 10), e.result.tm2 = c, e.result.len2 = h.len, e.result.gc2 = 100 * h.fgc, e.result.gc2 = Math.round(10 * e.result.gc2 /
            10);
          else {
            e.results2 = [], o.setCt(v).setMonosalt(r).setMethod(s).setDMSO(C);
            for (var N = 0; N < O.length; ++N) E = O[N], h = o.calcTm_nomm(E), c = Math.round(Math.round(10 * h.tm) / 10), e.results2.push({
              seq: E,
              tm: c,
              len: E.length,
              gc: Math.round(100 * h.fgc * 10 / 10)
            });
            c = e.results2.reduce((function(e, t) {
                return t.tm < e && (e = t.tm), e
              }), e.results2[0].tm), p = e.results2.reduce((function(e, t) {
                return t.tm > e && (e = t.tm), e
              }), e.results2[0].tm), e.result.tm2 = c !== p ? c + " - " + p : c, e.result.len2 = m.replace(/\s/g, "").length, e.result.gc2 = "---", e
              .result.tm2m = p
          }
        }
        if ("invalidseq" !== e.p1status && "invalidseq" !== e.p2status) {
          u = i.getAnnealTemp(g, l, m, c, b, $), e.result.ta = Math.round(u), console.log("checkpoint 3");
          var P = i.validateTm(l, c, u, $, b, t),
            q = i.validateBuffer($, b, t);
          if (x = [], w = [], P.hasCritWarnings && (x = P.critwarnings), P.hasWarnings && (w = P.warnings), q.hasCritWarnings && Array.prototype.push
            .apply(x, q.critwarnings), q.hasWarnings && Array.prototype.push.apply(w, q.warnings), e.result.itemlist = w, e.result.critlist = x, e
            .result && "---" === e.result.ta) e.showInfoMessage = !1;
          else switch (b) {
            case "Phusion":
            case "Phusion Hot Start Flex":
            case "Q5":
            case "Q5 Hot Start":
            case "Q5U Hot Start":
            case "Q5 Blood Direct":
              e.showInfoMessage = !0;
              break;
            case "Master Mix":
              0 === $.indexOf("phusion") || 0 === $.indexOf("q5") ? e.showInfoMessage = !0 : e.showInfoMessage = !1;
              break;
            default:
              e.showInfoMessage = !1
          }
          console.log("exiting runCalc2 normally")
        } else e.result.ta = "---"
      };
      e.$on("$routeChangeStart", (function() {
        ! function() {
          var n, r;
          for (n in console.log("saveState"), t.state || (t.state = {}, t.state.input = {}, t.state.result = {}, t.state.results1 = [], t.state
              .results2 = []), t.state.input || (t.state.input = {}), t.state.result || (t.state.result = {}), t.state.results1 || (t.state
              .results1 = []), t.state.results2 || (t.state.results2 = []), e.input) t.state.input[n] = e.input[n];
          for (n in e.result) t.state.result[n] = e.result[n];
          for (r in e.results1) t.state.results1[r] = e.results1[r];
          for (r in e.results2) t.state.results2[r] = e.results2[r];
          t.state.results1.showexpanded1 = e.results1.showexpanded1, t.state.results2.showexpanded2 = e.results2.showexpanded2, t.state.p1status =
            e.p1status = "", t.state.p2status = e.p2status = "", console.log(t.state)
        }()
      })), e.about = function() {
        u.path("#/about")
      }, e.switch2batch = function() {
        u.path("/batch")
      };
      ! function() {
        var n, r;
        for (n in console.log("loadState"), console.log(t.state), e.input || (e.input = {}), e.result || (e.result = {}), e.results1 || (e
          .results1 = []), e.results2 || (e.results2 = []), t.state.input) e.input[n] = t.state.input[n];
        for (n in t.state.result) e.result[n] = t.state.result[n];
        for (r in t.state.results1) e.results1[r] = t.state.results1[r];
        for (r in t.state.results2) e.results2[r] = t.state.results2[r];
        e.results1.showexpanded1 = t.state.results1.showexpanded1, e.results2.showexpanded2 = t.state.results2.showexpanded2, e.p1status = t.state
          .p1status, e.p2status = t.state.p2status, console.log(e.input), console.log(e.result)
      }(), 0 === e.input.groupKeys.length && (e.setGroups(), e.input.group = "Q5", e.input.productKeys = i.getProductKeysForGroup(e.input.group), e
          .input.products = i.getProductsForGroup(e.input.group), e.input.product = e.input.products[0].id, e.input.ct = 500), t.initScreenWidth <
        640 && t.screenMessage && (t.openmodal("screensize", "sm"), t.screenMessage = !1), e.showInfoMessage = !1, e.infoMessageText =
        "The annealing temperature for each polymerase is based on empirical observations of efficiency. The optimal annealing temperature for high fidelity hot start DNA polymerases like Q5 may differ significantly from that of Taq-based polymerases.",
        e.toggle_expanded2_display = function() {
          e.results2.showexpanded2 = !e.results2.showexpanded2, e.expanded2_toggle_label = e.results2.showexpanded2 ? "Hide" : "Show"
        }, e.toggle_expanded1_display = function() {
          e.results1.showexpanded1 = !e.results1.showexpanded1, e.expanded1_toggle_label = e.results1.showexpanded1 ? "Hide" : "Show"
        }, e.expanded2_toggle_label = e.results2.showexpanded2 ? "Hide" : "Show", e.expanded1_toggle_label = e.results1.showexpanded1 ? "Hide" : "Show"
    }
  }, {
    angular: 32
  }],
  9: [function(e, t, n) {
    "use strict"
  }, {}],
  10: [function(e, t, n) {
    "use strict";
    t.exports = function() {
      return {
        restrict: "A",
        replace: !0,
        scope: {
          processwith: "&"
        },
        template: '<div><div></div><input type="file" /></div>',
        link: function(e, t, n) {
          var r = angular.element(t.children().eq(1)),
            i = function(t) {
              var n = new FileReader;
              n.onload = function(t) {
                console.log("about to call processwith"), console.log(t.target.result), console.log(e.processwith), e.processwith({
                  filedata: t.target.result
                }), console.log("after call to processwith")
              }, n.readAsText((t.srcElement || t.target).files[0])
            };
          r.on("dragover", (function(e) {
            e.stopPropagation(), e.preventDefault(), e.dataTransfer.dropEffect = "copy"
          })), r.on("drop", i), r.on("change", i)
        }
      }
    }
  }, {}],
  11: [function(e, t, n) {
    "use strict";
    var r = e("angular").module("tmcApp");
    r.directive("fileInputTest", [e("./fileInputTest")]), r.directive("fileinput", [e("./fileinput")]), r.directive("onReadFile", ["$parse", e(
      "./onReadFile")])
  }, {
    "./fileInputTest": 9,
    "./fileinput": 10,
    "./onReadFile": 12,
    angular: 32
  }],
  12: [function(e, t, n) {
    "use strict";
    t.exports = function(e) {
      return {
        restrict: "A",
        scope: !1,
        link: function(t, n, r) {
          var i = e(r.onReadFile);
          n.on("change", (function(e) {
            var n = new FileReader;
            n.onload = function(e) {
              t.$apply((function() {
                i(t, {
                  filedata: e.target.result
                })
              }))
            }, n.readAsText((e.srcElement || e.target).files[0])
          }))
        }
      }
    }
  }, {}],
  13: [function(e, t, n) {
    "use strict";
    var r, i, o = {};
    (i = o).isNumeric = function(e) {
        return !isNaN(e / 1)
      }, i.trim = function(e) {
        for (var t = (e = e.replace(/^\s+/, "")).length - 1; t >= 0; t--)
          if (/\S/.test(e.charAt(t))) {
            e = e.slice(0, t + 1);
            break
          } return e
      }, i.DEGENERATES = {
        A: "A",
        C: "C",
        G: "G",
        T: "T",
        U: "U",
        W: "AT",
        S: "CG",
        M: "AC",
        K: "GT",
        R: "AG",
        Y: "CT",
        B: "CGT",
        D: "AGT",
        H: "ACT",
        V: "ACG",
        N: "ACGT",
        a: "a",
        c: "c",
        g: "g",
        t: "t",
        u: "u",
        w: "at",
        s: "cg",
        m: "ac",
        k: "gt",
        r: "ag",
        y: "ct",
        b: "cgt",
        d: "agt",
        h: "act",
        v: "acg",
        n: "acgt"
      }, i.COMPLEMENTS = function() {
        var e, t = "ACGTUNSWMKRYVDHBacgtunswmkryvdhb ",
          n = {};
        for (e = 0; e < t.length; ++e) n[t.charAt(e)] = "TGCAANSWKMYRBHDVtgcaanswkmyrbhdv ".charAt(e);
        return n
      }(), i.AMINO_ACID_TLC = {
        ALA: "A",
        ASX: "B",
        CYS: "C",
        ASP: "D",
        GLU: "E",
        PHE: "F",
        GLY: "G",
        HIS: "H",
        ILE: "I",
        LYS: "K",
        LEU: "L",
        MET: "M",
        ASN: "N",
        PYL: "O",
        PRO: "P",
        GLN: "Q",
        ARG: "R",
        SER: "S",
        THR: "T",
        SEC: "U",
        VAL: "V",
        TRP: "W",
        XAA: "X",
        TYR: "Y",
        GLX: "Z",
        TER: "*",
        Ala: "A",
        Asx: "B",
        Cys: "C",
        Asp: "D",
        Glu: "E",
        Phe: "F",
        Gly: "G",
        His: "H",
        Ile: "I",
        Lys: "K",
        Leu: "L",
        Met: "M",
        Asn: "N",
        Pyl: "O",
        Pro: "P",
        Gln: "Q",
        Arg: "R",
        Ser: "S",
        Thr: "T",
        Sec: "U",
        Val: "V",
        Trp: "W",
        Xaa: "X",
        Tyr: "Y",
        Glx: "Z",
        Ter: "*"
      }, i.AMINO_ACID_OLC = {
        A: "Ala",
        B: "Asx",
        C: "Cys",
        D: "Asp",
        E: "Glu",
        F: "Phe",
        G: "Gly",
        H: "His",
        I: "Ile",
        K: "Lys",
        L: "Leu",
        M: "Met",
        N: "Asn",
        O: "Pyl",
        P: "Pro",
        Q: "Gln",
        R: "Arg",
        S: "Ser",
        T: "Thr",
        U: "Sec",
        V: "Val",
        W: "Trp",
        X: "Xaa",
        Y: "Tyr",
        Z: "Glx",
        "*": "Ter"
      }, i.validAARegex = /([ABCDEFGHIKLMNOPQRSTUVWXYZ\*])/gi, i.invalidAARegex = /([^ABCDEFGHIKLMNOPQRSTUVWXYZ\*])/gi, i.validNARegex =
      /([ACGTUNSWMKRYVDHBacgtunswmkryvdhb])/g, i.invalidNARegex = /([^ACGTUNSWMKRYVDHBacgtunswmkryvdhb])/g, i.validNAStrictRegex = /([TCAGtcag])/g, i
      .invalidNAStrictRegex = /([^TCAGtcag])/g, i.isValidAASeq = function(e) {
        return !i.invalidAARegex.exec(e)
      }, i.isValidNASeq = function(e, t) {
        return !((t = t || !1) ? i.invalidNARegex : i.invalidNAStrictRegex).exec(e)
      }, i.aaOneToThree = function(e) {
        for (var t = e.split(""), n = [], r = 0; r < t.length; r++) {
          if (!i.AMINO_ACID_OLC[t[r]]) throw "Unknown one-letter-code: " + t[r];
          n.push(i.AMINO_ACID_OLC[t[r]])
        }
        return n.join("")
      }, i.CODONS = function() {
        var e, t, n, r, i, o, a = [],
          s = ["T", "C", "A", "G"];
        for (r = 0; r < 4; ++r)
          for (e = s[r], i = 0; i < 4; ++i)
            for (t = s[i], o = 0; o < 4; ++o) n = s[o], a.push(e + t + n);
        return a
      }(), i.CODON_TABLE_DATA = [
        "1\nFFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n---M---------------M---------------M----------------------------",
        "2\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG\n--------------------------------MMMM---------------M------------",
        "3\nFFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n----------------------------------MM----------------------------",
        "4\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n--MM---------------M------------MMMM---------------M------------",
        "5\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG\n---M----------------------------MMMM---------------M------------",
        "6\nFFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n-----------------------------------M----------------------------",
        "9\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\n-----------------------------------M---------------M------------",
        "10\nFFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n---M---------------M------------MMMM---------------M------------",
        "11\nFFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n---M---------------M------------MMMM---------------M------------",
        "12\nFFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n-------------------M---------------M----------------------------",
        "13\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG\n---M------------------------------MM---------------M------------",
        "14\nFFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\n-----------------------------------M----------------------------",
        "15\nFFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n-----------------------------------M----------------------------",
        "16\nFFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n-----------------------------------M----------------------------",
        "21\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG\n-----------------------------------M---------------M------------",
        "22\nFFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n-----------------------------------M----------------------------",
        "23\nFF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n--------------------------------M--M---------------M------------",
        "24\nFFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG\n---M---------------M---------------M---------------M------------",
        "25\nFFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n---M-------------------------------M---------------M------------",
        "90\nFFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\n---M---------------M------------MMMM---------------M------------"
      ], i.sanitizeSeq = function(e) {
        return e.replace(/[\s\d]+/gi, "")
      }, i.revcomp = function(e, t) {
        if ("" === (e = e || "")) return "";
        "undefined" !== t && null !== t || (t = !0);
        for (var n = e.split(""), r = 0; r < n.length; ++r)
          if (i.COMPLEMENTS[n[r]]) n[r] = i.COMPLEMENTS[n[r]];
          else if (t) throw "Unknown base: " + n[r];
        for (var o = "", a = n.length - 1; a >= 0; a--) o += n[a];
        return o
      }, i.rev = function(e) {
        if ("" === (e = e || "")) return "";
        for (var t = e.split(""), n = "", r = t.length - 1; r >= 0; r--) n += t[r];
        return n
      }, i.comp = function(e, t) {
        if ("" === (e = e || "")) return "";
        "undefined" !== t && null !== t || (t = !0);
        for (var n = e.split(""), r = 0; r < n.length; ++r)
          if (i.COMPLEMENTS[n[r]]) n[r] = i.COMPLEMENTS[n[r]];
          else if (t) throw "Unknown base: " + n[r];
        return n.join("")
      }, i.TranslationTable = function(e, t) {
        this.self = this, this.codons = e, this.starts = t, this.translate = function(n, r) {
          if (r = r || !1, 3 !== n.length) throw "Invalid codon: " + n;
          var o = n.split("").map((function(e) {
              return i.DEGENERATES[e].split("")
            })),
            a = function() {
              var e, t, n, r, i, a, s, u = [];
              for (i = 0; i < o[0].length; ++i) {
                for (r = [], e = o[0][i], a = 0; a < o[1].length; ++a)
                  for (t = o[1][a], s = 0; s < o[2].length; ++s) n = o[2][s], r.push(e + t + n);
                u.push(r)
              }
              return u
            }(),
            s = function() {
              for (var t = [], n = function(n) {
                  e[n] ? t.push(e[n]) : t.push("X")
                }, r = 0; r < a.length; r++) a[r].forEach(n);
              return t
            }(),
            u = function() {
              for (var e = [], n = function(n) {
                  t[n] ? e.push(t[n]) : e.push(!1)
                }, r = 0; r < a.length; r++) a[r].forEach(n);
              return e
            }();
          return [1 === s.length ? s.pop() : 2 === s.length ? s.indexOf("D") + s.indexOf("N") === 1 ? "B" : s.indexOf("E") + s.indexOf("Q") === 1 ?
            "Z" : s.indexOf("I") + s.indexOf("L") === 1 ? "J" : s[0] === s[1] ? s[0] : "X" : "X", !(u.length > 1) && u.pop()
          ]
        }
      }, i.buildhash = function(e, t) {
        var n = {};
        if (e.length !== t.length) throw "Arrays for buildhash must be equal length";
        for (var r = 0; r < e.length; r++) n[e[r].toString()] = t[r];
        return n
      }, i.CODON_TABLES = (r = {}, i.CODON_TABLE_DATA.forEach((function(e) {
        var t = e.split("\n");
        if (3 !== t.length) throw "codon table not parsed properly";
        var n = t[0],
          o = i.buildhash(i.CODONS, t[1].split("")),
          a = i.buildhash(i.CODONS, t[2].split("").map((function(e) {
            return "M" === e
          })));
        r[n] = new i.TranslationTable(o, a)
      })), r), i.translate = function(e, t, n, r, o) {
        if (n = n || "1", t = t || 0, r = r || [], o = o || [!0, !0], e = e.slice(t).toUpperCase().replace(/U/g, "T"), !i.CODON_TABLES[n])
        throw "Unknown translation table: " + n;
        if (e.length < 3) throw "Sequence length is too short.";
        var a, s = i.CODON_TABLES[n],
          u = s.translate(e.slice(0, 3), 3 === e.length),
          l = u[0],
          c = u[1];
        a = !o[0] && c ? ["M"] : o[0] || c || 0 !== t ? [l] : ["X"];
        for (var d = e.length, p = 3; p < d - d % 3; p += 3) a.push(s.translate(e.slice(p, p + 3), p === d - 3)[0]);
        return a = a.join(""), r.forEach((function(e) {
          if (e[0] % 3 != 0) throw "Invalid modification coordinate: " + e[0];
          var t = e[0] / 3;
          a = a.slice(0, t) + e[1] + a.slice(t)
        })), o[1] || "*" !== a[a.length - 1] || (a = a.slice(0, a.length - 1)), a
      }, i.terminalOverlapLen = function(e, t, n, r) {
        r = r || 0;
        var i, o = n = n || Math.min(e.length, t.length),
          a = function(t, n) {
            for (var r = 0, i = 0; i < e.length; ++i) t.charAt(i) !== n.charAt(i) && ++r;
            return r
          };
        for ((e.length < n || t.length < n) && (o = Math.min(e.length, t.length)), i = o; i > 0;) {
          if (a(e.slice(-i), t.slice(0, i)) <= r) return i;
          --i
        }
        return 0
      }, i.lcs = function(e, t) {
        if (!e || !t || 0 === e.length || 0 === t.length) return 0;
        var n, r, i = 0,
          o = 0,
          a = new Array(e.length);
        for (n = 0; n <= e.length; n++)
          for (a[n] = new Array(t.length), r = 0; r <= t.length; r++) a[n][r] = 0;
        for (var s = 0, u = e.length; s < u; s++)
          for (var l = 0, c = t.length; l < c; l++) e[s] === t[l] ? (a[s + 1][l + 1] = 0 === s || 0 === l ? 1 : a[s][l] + 1, a[s + 1][l + 1] > i && (i =
            a[s + 1][l + 1], o = s)) : a[s + 1][l + 1] = 0;
        return [i, e.slice(o - i + 1, o + 1)]
      }, i.LCS = function(e, t) {
        var n, r, i = e.length,
          o = t.length,
          a = [];
        for (n = 0; n <= i; n++) a.push([0]);
        for (r = 0; r < o; r++) a[0].push(0);
        for (n = 0; n < i; n++)
          for (r = 0; r < o; r++) a[n + 1][r + 1] = e[n] === t[r] ? a[n][r] + 1 : Math.max(a[n + 1][r], a[n][r + 1]);
        return function n(r, i) {
          return r * i == 0 ? "" : e[r - 1] === t[i - 1] ? n(r - 1, i - 1) + e[r - 1] : a[r][i - 1] > a[r - 1][i] ? n(r, i - 1) : n(r - 1, i)
        }(i, o)
      }, i.editDist = function(e, t) {
        if (0 === e.length) return t.length;
        if (0 === t.length) return e.length;
        var n, r, i = [];
        for (n = 0; n <= e.length; n++) i[n] = [n];
        for (r = 0; r <= t.length; r++) i[0][r] = r;
        for (n = 1; n <= e.length; n++)
          for (r = 1; r <= t.length; r++) e.charAt(n - 1) === t.charAt(r - 1) ? i[n][r] = i[n - 1][r - 1] : i[n][r] = Math.min(i[n - 1][r - 1] + 1, i[n][
            r - 1
          ] + 1, i[n - 1][r] + 1);
        return i[e.length][t.length]
      }, i.fGC = function(e) {
        for (var t = 0, n = 0; n < e.length; ++n) "g" !== e.charAt(n) && "G" !== e.charAt(n) && "C" !== e.charAt(n) && "c" !== e.charAt(n) || (t += 1);
        return t / e.length
      }, i.cntGC = function(e) {
        for (var t = 0, n = 0; n < e.length; ++n) "g" !== e.charAt(n) && "G" !== e.charAt(n) && "C" !== e.charAt(n) && "c" !== e.charAt(n) || (t += 1);
        return t
      }, i.randomDNASeq = function(e, t) {
        if (e < 0 && (e = 100), t && t.a && t.c && t.g && t.t && (t.a + t.c + t.g + t.t > 1 && (t.t = 1 - (t.a + t.c + t.g)), t.t < 0)) throw new Error(
          "Invalid base distribution");
        for (var n, r = ""; 0 < length;)(n = Math.random()) < t.a ? r += "a" : n < t.a + t.c ? r += "c" : n < t.a + t.c + t.g ? r += "g" : r += "t";
        return r
      }, i.replicateDNASeq = function(e, t) {
        for (var n = [], r = 0; r < t; ++r) n.push(e);
        return n
      }, i.cp = function(e) {
        for (var t, n, r = [], i = 0; i < e[0].length; ++i) t = e[0][i], r.push([t]);
        for (var o = 1; o < e.length;) {
          for (var a = [], s = 0; s < r.length; ++s) {
            n = r[s];
            for (var u = 0; u < e[o].length; ++u) t = e[o][u], a.push(n.concat([t]))
          }
          r = a.slice(), o++
        }
        return r
      }, i.expandDNASeq = function(e) {
        for (var t, n = /[^acgt]/gi, r = e.replace(/\s+/g, ""), o = []; null !== (t = n.exec(r));) o.push(t.index);
        for (var a, s = [], u = 0, l = 0; l < o.length; ++l)(a = o[l]) > u && s.push([r.slice(u, a)]), s.push(i.DEGENERATES[r[a]].split("")), u = a + 1;
        u < r.length && s.push([r.slice(u)]);
        for (var c, d = i.cp(s), p = [], f = 0; f < d.length; ++f) c = d[f], p.push(c.join(""));
        return p
      }, i.NASeq = function(e, t) {
        this._modelNASeq = {}, this._modelNASeq.stype = void 0, this._modelNASeq.isDoubleStranded = void 0, this._modelNASeq.topology = "linear", this
          ._modelNASeq.sid = void 0, this._modelNASeq.sname = void 0, this._modelNASeq.wseq = void 0, this._modelNASeq.cseq = void 0, this._modelNASeq
          .length = void 0, this._modelNASeq.sym = void 0, this._modelNASeq.gc_cnt = 0, this._modelNASeq.fgc = 0, this.setSeq(e), t && (t.stype && this
            .setSeqType(t.stype), t.topology && this.setSeqTopology(t.topology), t.sid && this.setSeqId(t.sid), t.sname && this.setSeqName(t.sname)),
          Object.defineProperty(this, "stype", {
            get: function() {
              return this._modelNASeq.stype
            },
            set: function(e) {
              this.setSeqType(e)
            }
          }), Object.defineProperty(this, "isDoubleStranded", {
            get: function() {
              return this._modelNASeq.isDoubleStranded
            },
            set: function(e) {
              this._modelNASeq.isDoubleStranded = !!e
            }
          }), Object.defineProperty(this, "sym", {
            get: function() {
              return this._modelNASeq.sym
            },
            set: function(e) {
              this._modelNASeq.sym = !!e
            }
          }), Object.defineProperty(this, "topology", {
            get: function() {
              return this._modelNASeq.topology
            },
            set: function(e) {
              this.setSeqTopology(e)
            }
          }), Object.defineProperty(this, "sid", {
            get: function() {
              return this._modelNASeq.sid
            },
            set: function(e) {
              this.setSeqId(e)
            }
          }), Object.defineProperty(this, "sname", {
            get: function() {
              return this._modelNASeq.sname
            },
            set: function(e) {
              this.setSeqName(e)
            }
          }), Object.defineProperty(this, "wseq", {
            get: function() {
              return this._modelNASeq.wseq
            },
            set: function(e) {
              this._modelNASeq.wseq = e
            }
          }), Object.defineProperty(this, "cseq", {
            get: function() {
              return this._modelNASeq.cseq
            },
            set: function(e) {
              this._modelNASeq.cseq = e
            }
          }), Object.defineProperty(this, "length", {
            get: function() {
              return this._modelNASeq.length
            },
            set: function(e) {
              this._modelNASeq.length = e
            }
          }), Object.defineProperty(this, "seq", {
            get: function() {
              return this._modelNASeq.wseq
            },
            set: function(e) {
              this.setSeq(e)
            }
          })
      }, i.NASeq.prototype = {
        setSeq: function(e) {
          if (this._modelNASeq.wseq = i.sanitizeSeq(e).toLowerCase(), !i.isValidNASeq(this._modelNASeq.wseq, !1)) throw new Error("Invalid seq");
          "" === this._modelNASeq.wseq ? this._modelNASeq.cseq = "" : this._modelNASeq.cseq = i.revcomp(this._modelNASeq.wseq), this._modelNASeq
            .length = this._modelNASeq.wseq.length, this._modelNASeq.sym = this._modelNASeq.wseq === this._modelNASeq.cseq, this._modelNASeq.gc_cnt = 0;
          for (var t = 0; t < this._modelNASeq.wseq.length; ++t) "g" !== this._modelNASeq.wseq.charAt(t) && "G" !== this._modelNASeq.wseq.charAt(t) &&
            "C" !== this._modelNASeq.wseq.charAt(t) && "c" !== this._modelNASeq.wseq.charAt(t) || (this._modelNASeq.gc_cnt += 1);
          return this._modelNASeq.length > 0 ? this._modelNASeq.fgc = this._modelNASeq.gc_cnt / this._modelNASeq.wseq.length : this._modelNASeq.fgc = 0,
            this
        },
        setSeqName: function(e) {
          return this._modelNASeq.sname = e, this
        },
        setSeqId: function(e) {
          return this._modelNASeq.sid = e, this
        },
        setSeqType: function(e) {
          return this._modelNASeq.stype = e, this
        },
        setSeqTopology: function(e) {
          return this._modelNASeq.topology = e, this
        },
        subSeq: function(e, t) {
          if (e < 0 && (e += this._modelNASeq.wseq.length), t < 0 && (t += this._modelNASeq.wseq.length), e < 0 || e >= this._modelNASeq.wseq.length ||
            t <= 0 || this._modelNASeq.wseq.length, e < t) return this._modelNASeq.wseq.slice(e, t);
          if ("circular" === this._modelNASeq.topology) return this._modelNASeq.wseq.slice(e) + this._modelNASeq.wseq.slice(0, t);
          throw new Error("startidx must be < endidx for linear seq")
        },
        splice: function(e, t, n) {
          if (n = n || "", arguments.length < 2) throw new Error("insufficient arguments");
          if (e < 0 && (e += this._modelNASeq.wseq.length), e < 0 || e >= this._modelNASeq.wseq.length) throw new Error("start index out of range");
          return e + t < this.wseq.length ? this.setSeq(this._modelNASeq.wseq.slice(0, e) + n + this._modelNASeq.wseq.slice(e + t)) : "linear" === this
            .topology ? this.setSeq(this._modelNASeq.wseq.slice(0, e) + n) : this.setSeq(this._modelNASeq.wseq.slice(e + t - this._modelNASeq.wseq
              .length, e) + n), this
        },
        replace: function(e, t, n) {
          if (n = n || "", arguments.length < 2) throw new Error("insufficient arguments");
          var r = t - e + 1;
          this.splice(e, r, n)
        },
        insert: function(e, t) {
          if (t < 0 || t >= this.wseq.length) throw new Error("start index out of range");
          this.setSeq(this._modelNASeq.wseq.slice(0, t + 1) + e + this._modelNASeq.wseq.slice(t + 1))
        },
        delete: function(e, t) {
          if (e < 0 || e > this._modelNASeq.wseq.length) throw new Error("start index out of range");
          if (t && (t < 0 || t > this._modelNASeq.wseq.length)) throw new Error("end index out of range");
          if (t)
            if (e > t) {
              if ("linear" === this.topology) throw new Error("end index smaller than start index");
              this.setSeq(this._modelNASeq.wseq.slice(t, e))
            } else this.setSeq(this._modelNASeq.wseq.slice(0, e) + this._modelNASeq.wseq.slice(t));
          else this.setSeq(this._modelNASeq.wseq.slice(0, e))
        },
        ligate: function(e) {
          return new i.NASeq(this._modelNASeq.wseq + e.wseq)
        },
        shift: function(e) {
          var t;
          e > 0 ? (t = this._modelNASeq.wseq.slice(e) + this._modelNASeq.wseq.slice(0, e), this.setSeq(t)) : e < 0 && (t = this._modelNASeq.wseq.slice(
            this.wseq.length + e) + this._modelNASeq.wseq.slice(0, this._modelNASeq.wseq.length + e), this.setSeq(t))
        },
        toJson: function() {
          return this._modelNASeq
        },
        expand: function() {
          var e = i.expandDNASeq(this.seq),
            t = this.stype,
            n = this.topology;
          return e.map((function(e) {
            new i.NASeq(e, {
              stype: t,
              topology: n
            })
          }))
        }
      }, void 0 !== n ? (void 0 !== t && t.exports && (n = t.exports = o), n.NEB_bioseq = o) : (void 0).NEB_bioseq = o, "function" == typeof define &&
      define.amd && define("NEB_bioseq", [], (function() {
        return o
      }))
  }, {}],
  14: [function(e, t, n) {
    "use strict";
    var r, i = {};
    (r = i).TmCalc = function(e) {
        this.init(e)
      }, r.TmCalc.prototype.nnBr = {
        aa: {
          dh: -9.1,
          ds: -22.2
        },
        tt: {
          dh: -9.1,
          ds: -22.2
        },
        at: {
          dh: -8.6,
          ds: -20.4
        },
        ta: {
          dh: -6,
          ds: -21.3
        },
        ca: {
          dh: -5.8,
          ds: -22.7
        },
        tg: {
          dh: -5.8,
          ds: -22.7
        },
        gt: {
          dh: -6.5,
          ds: -22.4
        },
        ac: {
          dh: -6.5,
          ds: -22.4
        },
        ct: {
          dh: -7.8,
          ds: -21
        },
        ag: {
          dh: -7.8,
          ds: -21
        },
        ga: {
          dh: -5.6,
          ds: -22.2
        },
        tc: {
          dh: -5.6,
          ds: -22.2
        },
        cg: {
          dh: -11.9,
          ds: -27.2
        },
        gc: {
          dh: -11.1,
          ds: -24.4
        },
        gg: {
          dh: -11,
          ds: -19.9
        },
        cc: {
          dh: -11,
          ds: -19.9
        }
      }, r.TmCalc.prototype.nnSa = {
        aa: {
          dh: -7.9,
          ds: -22.2
        },
        tt: {
          dh: -7.9,
          ds: -22.2
        },
        at: {
          dh: -7.2,
          ds: -20.4
        },
        ta: {
          dh: -7.2,
          ds: -21.3
        },
        ca: {
          dh: -8.5,
          ds: -22.7
        },
        tg: {
          dh: -8.5,
          ds: -22.7
        },
        gt: {
          dh: -8.4,
          ds: -22.4
        },
        ac: {
          dh: -8.4,
          ds: -22.4
        },
        ct: {
          dh: -7.8,
          ds: -21
        },
        ag: {
          dh: -7.8,
          ds: -21
        },
        ga: {
          dh: -8.2,
          ds: -22.2
        },
        tc: {
          dh: -8.2,
          ds: -22.2
        },
        cg: {
          dh: -10.6,
          ds: -27.2
        },
        gc: {
          dh: -9.8,
          ds: -24.4
        },
        gg: {
          dh: -8,
          ds: -19.9
        },
        cc: {
          dh: -8,
          ds: -19.9
        }
      }, r.TmCalc.prototype.nnde = {
        "-a/at": {
          dh: -.7,
          ds: -.8
        },
        "-a/ct": {
          dh: 4.4,
          ds: 14.9
        },
        "-a/gt": {
          dh: -1.6,
          ds: -3.6
        },
        "-a/tt": {
          dh: 2.9,
          ds: 10.4
        },
        "-c/ag": {
          dh: -2.1,
          ds: -3.9
        },
        "-c/cg": {
          dh: -.2,
          ds: -.1
        },
        "-c/gg": {
          dh: -3.9,
          ds: -11.2
        },
        "-c/tg": {
          dh: -4.4,
          ds: -13.1
        },
        "-g/ac": {
          dh: -5.9,
          ds: -16.5
        },
        "-g/cc": {
          dh: -2.6,
          ds: -7.4
        },
        "-g/gc": {
          dh: -3.2,
          ds: -10.4
        },
        "-g/tc": {
          dh: -5.2,
          ds: -15
        },
        "-t/aa": {
          dh: -.5,
          ds: -1.1
        },
        "-t/ca": {
          dh: 4.7,
          ds: 7.7
        },
        "-t/ga": {
          dh: -4.1,
          ds: -13.1
        },
        "-t/ta": {
          dh: -3.8,
          ds: -12.6
        },
        "a-/ta": {
          dh: -2.9,
          ds: -7.6
        },
        "a-/tc": {
          dh: -4.1,
          ds: -13
        },
        "a-/tg": {
          dh: -4.2,
          ds: -15
        },
        "a-/tt": {
          dh: -.2,
          ds: -.5
        },
        "aa/-t": {
          dh: .2,
          ds: 2.3
        },
        "aa/t-": {
          dh: -.5,
          ds: -1.1
        },
        "ac/-g": {
          dh: -6.3,
          ds: -17.1
        },
        "ac/t-": {
          dh: 4.7,
          ds: 7.7
        },
        "ag/-c": {
          dh: -3.7,
          ds: -10
        },
        "ag/t-": {
          dh: -4.1,
          ds: -13.1
        },
        "at/-a": {
          dh: -2.9,
          ds: -7.6
        },
        "at/t-": {
          dh: -3.8,
          ds: -12.6
        },
        "c-/ga": {
          dh: -3.7,
          ds: -10
        },
        "c-/gc": {
          dh: -4,
          ds: -11.9
        },
        "c-/gg": {
          dh: -3.9,
          ds: -10.9
        },
        "c-/gt": {
          dh: -4.9,
          ds: -13.8
        },
        "ca/-t": {
          dh: .6,
          ds: 3.3
        },
        "ca/g-": {
          dh: -5.9,
          ds: -16.5
        },
        "cc/-g": {
          dh: -4.4,
          ds: -12.6
        },
        "cc/g-": {
          dh: -2.6,
          ds: -7.4
        },
        "cg/-c": {
          dh: -4,
          ds: -11.9
        },
        "cg/g-": {
          dh: -3.2,
          ds: -10.4
        },
        "ct/-a": {
          dh: -4.1,
          ds: -13
        },
        "ct/g-": {
          dh: -5.2,
          ds: -15
        },
        "g-/ca": {
          dh: -6.3,
          ds: -17.1
        },
        "g-/cc": {
          dh: -4.4,
          ds: -12.6
        },
        "g-/cg": {
          dh: -5.1,
          ds: -14
        },
        "g-/ct": {
          dh: -4,
          ds: -10.9
        },
        "ga/-t": {
          dh: -1.1,
          ds: -1.6
        },
        "ga/c-": {
          dh: -2.1,
          ds: -3.9
        },
        "gc/-g": {
          dh: -5.1,
          ds: -14
        },
        "gc/c-": {
          dh: -.2,
          ds: -.1
        },
        "gg/-c": {
          dh: -3.9,
          ds: -10.9
        },
        "gg/c-": {
          dh: -3.9,
          ds: -11.2
        },
        "gt/-a": {
          dh: -4.2,
          ds: -15
        },
        "gt/c-": {
          dh: -4.4,
          ds: -13.1
        },
        "t-/aa": {
          dh: .2,
          ds: 2.3
        },
        "t-/ac": {
          dh: .6,
          ds: 3.3
        },
        "t-/ag": {
          dh: -1.1,
          ds: -1.6
        },
        "t-/at": {
          dh: -6.9,
          ds: -20
        },
        "ta/-t": {
          dh: -6.9,
          ds: -20
        },
        "ta/a-": {
          dh: -.7,
          ds: -.8
        },
        "tc/-g": {
          dh: -4,
          ds: -10.9
        },
        "tc/a-": {
          dh: 4.4,
          ds: 14.9
        },
        "tg/-c": {
          dh: -4.9,
          ds: -13.8
        },
        "tg/a-": {
          dh: -1.6,
          ds: -3.6
        },
        "tt/-a": {
          dh: -.2,
          ds: -.5
        },
        "tt/a-": {
          dh: 2.9,
          ds: 10.4
        }
      }, r.TmCalc.prototype.nndhds = {
        "aa/at": {
          dh: 4.7,
          ds: 12.9
        },
        "aa/ct": {
          dh: 7.6,
          ds: 20.2
        },
        "aa/gt": {
          dh: 3,
          ds: 7.4
        },
        "aa/ta": {
          dh: 1.2,
          ds: 1.7
        },
        "aa/tc": {
          dh: -9.1,
          ds: -24
        },
        "aa/tg": {
          dh: -.6,
          ds: -2.3
        },
        "aa/tt": {
          dh: -7.6,
          ds: -21.3
        },
        "ac/ag": {
          dh: -2.9,
          ds: -9.8
        },
        "ac/cg": {
          dh: -.7,
          ds: -3.8
        },
        "ac/gg": {
          dh: .5,
          ds: 3.2
        },
        "ac/ta": {
          dh: 5.3,
          ds: 14.6
        },
        "ac/tc": {
          dh: 0,
          ds: -4.4
        },
        "ac/tg": {
          dh: -8.4,
          ds: -22.4
        },
        "ac/tt": {
          dh: .7,
          ds: .2
        },
        "ag/ac": {
          dh: -.9,
          ds: -4.2
        },
        "ag/cc": {
          dh: .6,
          ds: -.6
        },
        "ag/gc": {
          dh: -4,
          ds: -13.2
        },
        "ag/ta": {
          dh: -.7,
          ds: -2.3
        },
        "ag/tc": {
          dh: -7.8,
          ds: -21
        },
        "ag/tg": {
          dh: -3.1,
          ds: -9.5
        },
        "ag/tt": {
          dh: 1,
          ds: .9
        },
        "at/aa": {
          dh: 1.2,
          ds: 1.7
        },
        "at/ca": {
          dh: 5.3,
          ds: 14.6
        },
        "at/ga": {
          dh: -.7,
          ds: -2.3
        },
        "at/ta": {
          dh: -7.2,
          ds: -20.4
        },
        "at/tc": {
          dh: -1.2,
          ds: -6.2
        },
        "at/tg": {
          dh: -2.5,
          ds: -8.3
        },
        "at/tt": {
          dh: -2.7,
          ds: -10.8
        },
        "ca/at": {
          dh: 3.4,
          ds: 8
        },
        "ca/ct": {
          dh: 6.1,
          ds: 16.4
        },
        "ca/ga": {
          dh: -.9,
          ds: -4.2
        },
        "ca/gc": {
          dh: 1.9,
          ds: 3.7
        },
        "ca/gg": {
          dh: -.7,
          ds: -2.3
        },
        "ca/gt": {
          dh: -8.5,
          ds: -22.7
        },
        "ca/tt": {
          dh: 1,
          ds: .7
        },
        "cc/ag": {
          dh: 5.2,
          ds: 14.2
        },
        "cc/cg": {
          dh: 3.6,
          ds: 8.9
        },
        "cc/ga": {
          dh: .6,
          ds: -.6
        },
        "cc/gc": {
          dh: -1.5,
          ds: -7.2
        },
        "cc/gg": {
          dh: -8,
          ds: -19.9
        },
        "cc/gt": {
          dh: -.8,
          ds: -4.5
        },
        "cc/tg": {
          dh: 5.2,
          ds: 13.5
        },
        "cg/ac": {
          dh: 1.9,
          ds: 3.7
        },
        "cg/cc": {
          dh: -1.5,
          ds: -7.2
        },
        "cg/ga": {
          dh: -4,
          ds: -13.2
        },
        "cg/gc": {
          dh: -10.6,
          ds: -27.2
        },
        "cg/gg": {
          dh: -4.9,
          ds: -15.3
        },
        "cg/gt": {
          dh: -4.1,
          ds: -11.7
        },
        "cg/tc": {
          dh: -1.5,
          ds: -6.1
        },
        "ct/aa": {
          dh: -9.1,
          ds: -24
        },
        "ct/ca": {
          dh: 0,
          ds: -4.4
        },
        "ct/ga": {
          dh: -7.8,
          ds: -21
        },
        "ct/gc": {
          dh: -1.5,
          ds: -6.1
        },
        "ct/gg": {
          dh: -2.8,
          ds: -8
        },
        "ct/gt": {
          dh: -5,
          ds: -15.8
        },
        "ct/ta": {
          dh: -1.2,
          ds: -6.2
        },
        "ga/at": {
          dh: .7,
          ds: .7
        },
        "ga/ca": {
          dh: -2.9,
          ds: -9.8
        },
        "ga/cc": {
          dh: 5.2,
          ds: 14.2
        },
        "ga/cg": {
          dh: -.6,
          ds: -1
        },
        "ga/ct": {
          dh: -8.2,
          ds: -22.2
        },
        "ga/gt": {
          dh: 1.6,
          ds: 3.6
        },
        "ga/tt": {
          dh: -1.3,
          ds: -5.3
        },
        "gc/ag": {
          dh: -.6,
          ds: -1
        },
        "gc/ca": {
          dh: -.7,
          ds: -3.8
        },
        "gc/cc": {
          dh: 3.6,
          ds: 8.9
        },
        "gc/cg": {
          dh: -9.8,
          ds: -24.4
        },
        "gc/ct": {
          dh: 2.3,
          ds: 5.4
        },
        "gc/gg": {
          dh: -6,
          ds: -15.8
        },
        "gc/tg": {
          dh: -4.4,
          ds: -12.3
        },
        "gg/ac": {
          dh: -.7,
          ds: -2.3
        },
        "gg/ca": {
          dh: .5,
          ds: 3.2
        },
        "gg/cc": {
          dh: -8,
          ds: -19.9
        },
        "gg/cg": {
          dh: -6,
          ds: -15.8
        },
        "gg/ct": {
          dh: 3.3,
          ds: 10.4
        },
        "gg/gc": {
          dh: -4.9,
          ds: -15.3
        },
        "gg/tc": {
          dh: -2.8,
          ds: -8
        },
        "gt/aa": {
          dh: -.6,
          ds: -2.3
        },
        "gt/ca": {
          dh: -8.4,
          ds: -22.4
        },
        "gt/cc": {
          dh: 5.2,
          ds: 13.5
        },
        "gt/cg": {
          dh: -4.4,
          ds: -12.3
        },
        "gt/ct": {
          dh: -2.2,
          ds: -8.4
        },
        "gt/ga": {
          dh: -3.1,
          ds: -9.5
        },
        "gt/ta": {
          dh: -2.5,
          ds: -8.3
        },
        "ta/aa": {
          dh: 4.7,
          ds: 12.9
        },
        "ta/ac": {
          dh: 3.4,
          ds: 8
        },
        "ta/ag": {
          dh: .7,
          ds: .7
        },
        "ta/at": {
          dh: -7.2,
          ds: -20.4
        },
        "ta/ct": {
          dh: 1.2,
          ds: .7
        },
        "ta/gt": {
          dh: -.1,
          ds: -1.7
        },
        "ta/tt": {
          dh: .2,
          ds: -1.5
        },
        "tc/aa": {
          dh: 7.6,
          ds: 20.2
        },
        "tc/ac": {
          dh: 6.1,
          ds: 16.4
        },
        "tc/ag": {
          dh: -8.2,
          ds: -22.2
        },
        "tc/at": {
          dh: 1.2,
          ds: .7
        },
        "tc/cg": {
          dh: 2.3,
          ds: 5.4
        },
        "tc/gg": {
          dh: 3.3,
          ds: 10.4
        },
        "tc/tg": {
          dh: -2.2,
          ds: -8.4
        },
        "tg/aa": {
          dh: 3,
          ds: 7.4
        },
        "tg/ac": {
          dh: -8.5,
          ds: -22.7
        },
        "tg/ag": {
          dh: 1.6,
          ds: 3.6
        },
        "tg/at": {
          dh: -.1,
          ds: -1.7
        },
        "tg/cc": {
          dh: -.8,
          ds: -4.5
        },
        "tg/gc": {
          dh: -4.1,
          ds: -11.7
        },
        "tg/tc": {
          dh: -5,
          ds: -15.8
        },
        "tt/aa": {
          dh: -7.6,
          ds: -21.3
        },
        "tt/ac": {
          dh: 1,
          ds: .7
        },
        "tt/ag": {
          dh: -1.3,
          ds: -5.3
        },
        "tt/at": {
          dh: .2,
          ds: -1.5
        },
        "tt/ca": {
          dh: .7,
          ds: .2
        },
        "tt/ga": {
          dh: 1,
          ds: .9
        },
        "tt/ta": {
          dh: -2.7,
          ds: -10.8
        },
        "ai/tc": {
          dh: -8.9,
          ds: -25.5
        },
        "ti/ac": {
          dh: -5.9,
          ds: -17.4
        },
        "ac/ti": {
          dh: -8.8,
          ds: -25.4
        },
        "tc/ai": {
          dh: -4.9,
          ds: -13.9
        },
        "ci/gc": {
          dh: -5.4,
          ds: -13.7
        },
        "gi/cc": {
          dh: -6.8,
          ds: -19.1
        },
        "cc/gi": {
          dh: -8.3,
          ds: -23.8
        },
        "gc/ci": {
          dh: -5,
          ds: -12.6
        },
        "ai/ta": {
          dh: -8.3,
          ds: -25
        },
        "ti/aa": {
          dh: -3.4,
          ds: -11.2
        },
        "aa/ti": {
          dh: -.7,
          ds: -2.6
        },
        "ta/ai": {
          dh: -1.3,
          ds: -4.6
        },
        "ci/ga": {
          dh: 2.6,
          ds: 8.9
        },
        "gi/ca": {
          dh: -7.8,
          ds: -21.1
        },
        "ca/gi": {
          dh: -7,
          ds: -20
        },
        "ga/ci": {
          dh: -7.6,
          ds: -20.2
        },
        "ai/tt": {
          dh: .49,
          ds: -.7
        },
        "ti/at": {
          dh: -6.5,
          ds: -22
        },
        "at/ti": {
          dh: -5.6,
          ds: -18.7
        },
        "tt/ai": {
          dh: -.8,
          ds: -4.3
        },
        "ci/gt": {
          dh: -1,
          ds: -2.4
        },
        "gi/ct": {
          dh: -3.5,
          ds: -10.6
        },
        "ct/gi": {
          dh: .1,
          ds: -1
        },
        "gt/ci": {
          dh: -4.3,
          ds: -12.1
        },
        "ai/tg": {
          dh: -4.9,
          ds: -15.8
        },
        "ti/ag": {
          dh: -1.9,
          ds: -8.5
        },
        "ag/ti": {
          dh: .1,
          ds: -1.8
        },
        "tg/ai": {
          dh: 1,
          ds: 1
        },
        "ci/gg": {
          dh: 7.1,
          ds: 21.3
        },
        "gi/cg": {
          dh: -1.1,
          ds: -3.2
        },
        "cg/gi": {
          dh: 5.8,
          ds: 16.9
        },
        "gg/ci": {
          dh: -7.6,
          ds: -22
        },
        "ai/ti": {
          dh: -3.3,
          ds: -11.9
        },
        "ti/ai": {
          dh: .1,
          ds: -2.3
        },
        "ci/gi": {
          dh: 1.3,
          ds: 3
        },
        "gi/ci": {
          dh: -.5,
          ds: -1.3
        }
      }, r.TmCalc.prototype.nnloop = {
        3: {
          dh: 0,
          ds: -10.3
        },
        4: {
          dh: 0,
          ds: -11.6
        },
        5: {
          dh: 0,
          ds: -12.9
        },
        6: {
          dh: 0,
          ds: -14.2
        },
        7: {
          dh: 0,
          ds: -14.8
        },
        8: {
          dh: 0,
          ds: -15.5
        },
        9: {
          dh: 0,
          ds: -15.8
        },
        10: {
          dh: 0,
          ds: -15.8
        },
        12: {
          dh: 0,
          ds: -16.8
        },
        14: {
          dh: 0,
          ds: -17.4
        },
        16: {
          dh: 0,
          ds: -18.1
        },
        18: {
          dh: 0,
          ds: -18.7
        },
        20: {
          dh: 0,
          ds: -19
        },
        25: {
          dh: 0,
          ds: -20.3
        },
        30: {
          dh: 0,
          ds: -21.3
        }
      }, r.TmCalc.prototype.nnbulge = {
        1: {
          dh: 0,
          ds: -12.9
        },
        2: {
          dh: 0,
          ds: -9.4
        },
        3: {
          dh: 0,
          ds: -10
        },
        4: {
          dh: 0,
          ds: -10.3
        },
        5: {
          dh: 0,
          ds: -10.6
        },
        6: {
          dh: 0,
          ds: -11.3
        },
        7: {
          dh: 0,
          ds: -11.9
        },
        8: {
          dh: 0,
          ds: -12.6
        },
        9: {
          dh: 0,
          ds: -13.2
        },
        10: {
          dh: 0,
          ds: -13.9
        },
        12: {
          dh: 0,
          ds: -14.5
        },
        14: {
          dh: 0,
          ds: -15.5
        },
        16: {
          dh: 0,
          ds: -16.1
        },
        18: {
          dh: 0,
          ds: -16.8
        },
        20: {
          dh: 0,
          ds: -17.1
        },
        25: {
          dh: 0,
          ds: -28.1
        },
        30: {
          dh: 0,
          ds: -19
        }
      }, r.TmCalc.prototype.nntmm = {
        "aa/ta": {
          dh: -3.1,
          ds: -7.8
        },
        "aa/tc": {
          dh: -1.6,
          ds: -4
        },
        "aa/tg": {
          dh: -1.9,
          ds: -4.4
        },
        "ac/ta": {
          dh: -1.8,
          ds: -3.8
        },
        "ac/tc": {
          dh: -.1,
          ds: .5
        },
        "ac/tt": {
          dh: -.9,
          ds: -1.7
        },
        "ag/ta": {
          dh: -2.5,
          ds: -5.9
        },
        "ag/tg": {
          dh: -1.1,
          ds: -2.1
        },
        "ag/tt": {
          dh: -3.2,
          ds: -8.7
        },
        "at/tc": {
          dh: -2.3,
          ds: -6.3
        },
        "at/tg": {
          dh: -3.5,
          ds: -9.4
        },
        "at/tt": {
          dh: -2.4,
          ds: -6.5
        },
        "ca/ga": {
          dh: -4.3,
          ds: -10.7
        },
        "ca/gc": {
          dh: -2.6,
          ds: -5.9
        },
        "ca/gg": {
          dh: -3.9,
          ds: -9.6
        },
        "cc/ga": {
          dh: -2.7,
          ds: -6
        },
        "cc/gc": {
          dh: -2.1,
          ds: -5.1
        },
        "cc/gt": {
          dh: -3.2,
          ds: -8
        },
        "cg/ga": {
          dh: -6,
          ds: -15.5
        },
        "cg/gg": {
          dh: -3.8,
          ds: -9.5
        },
        "cg/gt": {
          dh: -3.8,
          ds: -9
        },
        "ct/gc": {
          dh: -3.9,
          ds: -10.6
        },
        "ct/gg": {
          dh: -6.6,
          ds: -18.7
        },
        "ct/gt": {
          dh: -6.1,
          ds: -16.9
        },
        "ga/ca": {
          dh: -8,
          ds: -22.5
        },
        "ga/cc": {
          dh: -5,
          ds: -13.8
        },
        "ga/cg": {
          dh: -4.3,
          ds: -11.1
        },
        "gc/ca": {
          dh: -3.2,
          ds: -7.1
        },
        "gc/cc": {
          dh: -3.9,
          ds: -10.6
        },
        "gc/ct": {
          dh: -4.9,
          ds: -13.5
        },
        "gg/ca": {
          dh: -4.6,
          ds: -11.4
        },
        "gg/cg": {
          dh: -.7,
          ds: -19.2
        },
        "gg/ct": {
          dh: -5.7,
          ds: -15.9
        },
        "gt/cc": {
          dh: -3,
          ds: -7.8
        },
        "gt/cg": {
          dh: -5.9,
          ds: -16.1
        },
        "gt/ct": {
          dh: -7.4,
          ds: -21.2
        },
        "ta/aa": {
          dh: -2.5,
          ds: -6.3
        },
        "ta/ac": {
          dh: -2.3,
          ds: -5.9
        },
        "ta/ag": {
          dh: -2,
          ds: -4.7
        },
        "tc/aa": {
          dh: -2.7,
          ds: -7
        },
        "tc/ac": {
          dh: -.7,
          ds: -1.3
        },
        "tc/at": {
          dh: -2.5,
          ds: -6.3
        },
        "tg/aa": {
          dh: -2.4,
          ds: -5.8
        },
        "tg/ag": {
          dh: -1.1,
          ds: -2.7
        },
        "tg/at": {
          dh: -3.9,
          ds: -10.5
        },
        "tt/ac": {
          dh: -.7,
          ds: -1.2
        },
        "tt/ag": {
          dh: -3.6,
          ds: -9.8
        },
        "tt/at": {
          dh: -3.2,
          ds: -8.9
        },
        "at/aa": {
          dh: -3.1,
          ds: -7.8
        },
        "ct/aa": {
          dh: -1.6,
          ds: -4
        },
        "gt/aa": {
          dh: -1.9,
          ds: -4.4
        },
        "at/ca": {
          dh: -1.8,
          ds: -3.8
        },
        "ct/ca": {
          dh: -.1,
          ds: .5
        },
        "tt/ca": {
          dh: -.9,
          ds: -1.7
        },
        "at/ga": {
          dh: -2.5,
          ds: -5.9
        },
        "gt/ga": {
          dh: -1.1,
          ds: -2.1
        },
        "tt/ga": {
          dh: -3.2,
          ds: -8.7
        },
        "ct/ta": {
          dh: -2.3,
          ds: -6.3
        },
        "gt/ta": {
          dh: -3.5,
          ds: -9.4
        },
        "tt/ta": {
          dh: -2.4,
          ds: -6.5
        },
        "ag/ac": {
          dh: -4.3,
          ds: -10.7
        },
        "cg/ac": {
          dh: -2.6,
          ds: -5.9
        },
        "gg/ac": {
          dh: -3.9,
          ds: -9.6
        },
        "ag/cc": {
          dh: -2.7,
          ds: -6
        },
        "cg/cc": {
          dh: -2.1,
          ds: -5.1
        },
        "tg/cc": {
          dh: -3.2,
          ds: -8
        },
        "ag/gc": {
          dh: -6,
          ds: -15.5
        },
        "gg/gc": {
          dh: -3.8,
          ds: -9.5
        },
        "tg/gc": {
          dh: -3.8,
          ds: -9
        },
        "cg/tc": {
          dh: -3.9,
          ds: -10.6
        },
        "gg/tc": {
          dh: -6.6,
          ds: -18.7
        },
        "tg/tc": {
          dh: -6.1,
          ds: -16.9
        },
        "ac/ag": {
          dh: -8,
          ds: -22.5
        },
        "cc/ag": {
          dh: -5,
          ds: -13.8
        },
        "gc/ag": {
          dh: -4.3,
          ds: -11.1
        },
        "ac/cg": {
          dh: -3.2,
          ds: -7.1
        },
        "cc/cg": {
          dh: -3.9,
          ds: -10.6
        },
        "tc/cg": {
          dh: -4.9,
          ds: -13.5
        },
        "ac/gg": {
          dh: -4.6,
          ds: -11.4
        },
        "gc/gg": {
          dh: -.7,
          ds: -19.2
        },
        "tc/gg": {
          dh: -5.7,
          ds: -15.9
        },
        "cc/tg": {
          dh: -3,
          ds: -7.8
        },
        "gc/tg": {
          dh: -5.9,
          ds: -16.1
        },
        "tc/tg": {
          dh: -7.4,
          ds: -21.2
        },
        "aa/at": {
          dh: -2.5,
          ds: -6.3
        },
        "ca/at": {
          dh: -2.3,
          ds: -5.9
        },
        "ga/at": {
          dh: -2,
          ds: -4.7
        },
        "aa/ct": {
          dh: -2.7,
          ds: -7
        },
        "ca/ct": {
          dh: -.7,
          ds: -1.3
        },
        "ta/ct": {
          dh: -2.5,
          ds: -6.3
        },
        "aa/gt": {
          dh: -2.4,
          ds: -5.8
        },
        "ga/gt": {
          dh: -1.1,
          ds: -2.7
        },
        "ta/gt": {
          dh: -3.9,
          ds: -10.5
        },
        "ca/tt": {
          dh: -.7,
          ds: -1.2
        },
        "ga/tt": {
          dh: -3.6,
          ds: -9.8
        },
        "ta/tt": {
          dh: -3.2,
          ds: -8.9
        }
      }, r.TmCalc.prototype.R = 1.987, r.TmCalc.prototype.dSBr = {
        aa: -24,
        tt: -24,
        at: -23.9,
        ta: -16.9,
        ca: -12.9,
        tg: -12.9,
        gt: -17.3,
        ac: -17.3,
        ct: -20.8,
        ag: -20.8,
        ga: -13.5,
        tc: -13.5,
        cg: -27.8,
        gc: -26.7,
        gg: -26.6,
        cc: -26.6
      }, r.TmCalc.prototype.dSSa = {
        aa: -22.2,
        tt: -22.2,
        at: -20.4,
        ta: -21.3,
        ca: -22.7,
        tg: -22.7,
        gt: -22.4,
        ac: -22.4,
        ct: -21,
        ag: -21,
        ga: -22.2,
        tc: -22.2,
        cg: -27.2,
        gc: -24.4,
        gg: -19.9,
        cc: -19.9
      }, r.TmCalc.prototype.dHBr = {
        aa: -9.1,
        tt: -9.1,
        at: -8.6,
        ta: -6,
        ca: -5.8,
        tg: -5.8,
        gt: -6.5,
        ac: -6.5,
        ct: -7.8,
        ag: -7.8,
        ga: -5.6,
        tc: -5.6,
        cg: -11.9,
        gc: -11.1,
        gg: -11,
        cc: -11
      }, r.TmCalc.prototype.dHSa = {
        aa: -7.9,
        tt: -7.9,
        at: -7.2,
        ta: -7.2,
        ca: -8.5,
        tg: -8.5,
        gt: -8.4,
        ac: -8.4,
        ct: -7.8,
        ag: -7.8,
        ga: -8.2,
        tc: -8.2,
        cg: -10.6,
        gc: -9.8,
        gg: -8,
        cc: -8
      }, r.TmCalc.prototype.init = function(e) {
        e ? e.seq && e.compseq ? this.setSeq(e.seq, e.compseq) : e.seq ? this.setSeq(e.seq) : (this.wseq = "", this.cseq = "") : (this.wseq = "", this
            .cseq = "", e = {}), this.setCt(e.ct || .25), this.setMonosalt(e.monosalt || 50), this.setDisalt(e.disalt || 0), this.setDMSO(e.dmso || 0),
          this.setMethod(e.method || 4), this.pairmap = e.pairmap || ""
      }, r.TmCalc.prototype.saltCorrect = function() {
        var e, t, n, r;
        if (0 === this.saltc && 0 === this.disaltc) throw {
          name: "Missing salt",
          message: "Unable to calculate salt correction because mono and divalent salt are missing."
        };
        this.sc_sch = 16.6 * Math.log(this.saltc) / Math.LN10, this.sc_sl = .368 * this.wseq.length * Math.log(this.saltc), t = 3.92, n = 1.42, r = 8.31,
          (e = this.saltc > 0 ? Math.sqrt(this.disaltc) / this.saltc : -1) >= 0 && e < .22 ? this.sc_ow = 1e-5 * (4.29 * this.fgc - 3.95) * Math.log(this
            .saltc) + 94e-7 * Math.log(this.saltc) * Math.log(this.saltc) : (e >= .22 && e < 6 && (t = 3.92 * (.843 - .352 * Math.sqrt(this.saltc) * Math
            .log(this.saltc)), n = 1.42 * (1.279 - .00403 * Math.log(this.saltc) - .00803 * Math.log(this.saltc) * Math.log(this.saltc)), r = 8.31 * (
            .486 - .258 * Math.log(this.saltc) + .00525 * Math.log(this.saltc) * Math.log(this.saltc) * Math.log(this.saltc))), this.sc_ow = 1e-5 * (t +
            -.911 * Math.log(this.disaltc) + this.fgc * (6.26 + n * Math.log(this.disaltc)) + .5 / (this.wseq.length - 1) * (52.5 * Math.log(this
              .disaltc) - 48.2 + r * Math.log(this.disaltc) * Math.log(this.disaltc))))
      }, r.TmCalc.prototype.setMethod = function(e) {
        return this.method = e, this
      }, r.TmCalc.prototype.setSeq = function(e, t) {
        var n;
        for (this.wseq = e.toLowerCase().replace(/\s/g, "").replace(/u/g, "t"), t ? this.cseq = t.toLowerCase().replace(/\s/g, "") : (this.cseq = this
            .revcomp(this.wseq, !1), this.pairmap = Array(this.wseq.length + 1).join("|")), this.sym = this.wseq === this.cseq, this.gc_cnt = 0, this
          .at_cnt = 0, n = 0; n < this.wseq.length; ++n) "g" === this.wseq.charAt(n) || "c" === this.wseq.charAt(n) ? this.gc_cnt += 1 : "a" !== this.wseq
          .charAt(n) && "t" !== this.wseq.charAt(n) || (this.at_cnt += 1);
        return this.fgc = this.gc_cnt / (this.at_cnt + this.gc_cnt), this.saltCorrect(), this
      }, r.TmCalc.prototype.setCt = function(e) {
        return this.ct = e, this.primerc = 1e-6 * this.ct, this
      }, r.TmCalc.prototype.setMonosalt = function(e) {
        return this.salt = e, this.saltc = .001 * this.salt, this.saltCorrect(), this
      }, r.TmCalc.prototype.setDisalt = function(e) {
        return this.disalt = e, this.disaltc = .001 * this.disalt, this.saltCorrect(), this
      }, r.TmCalc.prototype.setDMSO = function(e) {
        return this.dmso = e, this
      }, r.TmCalc.prototype.buildPairMap = function() {
        var e, t, n = this.wseq.split(""),
          r = this.cseq.split("").reverse(),
          i = function() {
            var e, t = "ACGTUNSWMKRYVDHBacgtunswmkryvdhb ",
              n = {};
            for (e = 0; e < t.length; ++e) n[t.charAt(e)] = "TGCAANSWKMYRBHDVtgcaanswkmyrbhdv ".charAt(e);
            return n
          }();
        for (e = 0; e < n.length; ++e) n[e] === i(r[e]) ? this.pairmap += "|" : "-" === n[e] || "-" === r[e] ? this.pairmap += "^" : this.pairmap += "x";
        for (var o = /x{3,}/g; null !== o.exec(this.pairmap););
        for (t = this.pairmap.split(""), e = 0; e < t.length - 2; ++e) "x" !== t[e] && "_" !== t[e] || "x" !== t[e + 1] && "_" !== t[e + 1] || "x" !== t[
          e + 2] && "_" !== t[e + 2] || (t[e] = "_", t[e + 1] = "_", t[e + 2] = "_");
        for (this.pairmap = t.join(""); - 1 !== this.pairmap.search(/\^+_+|_+\^+/);) this.pairmap = this.pairmap.replace(/\^+_+|_+\^+/, (function(e, t,
        n) {
          n.slice(0, t);
          var r = new Array(e.length + 1).join("_");
          return n.slice(t + e.length), r
        }));
        for (; - 1 !== this.pairmap.search(/\^+x|x\^+/);) this.pairmap = this.pairmap.replace(/\^+x|x\^+/, (function(e, t, n) {
          n.slice(0, t);
          var r = new Array(e.length + 1).join("_");
          return n.slice(t + e.length), r
        }));
        for (; - 1 !== this.pairmap.search(/_+[x^]|[x^]_+/);) this.pairmap = this.pairmap.replace(/_+[x^]|[x^]_+/, (function(e, t, n) {
          n.slice(0, t);
          var r = new Array(e.length + 1).join("_");
          return n.slice(t + e.length), r
        }))
      }, r.TmCalc.prototype.revcomp = function e(t) {
        if (e.COMPLEMENTS || (e.COMPLEMENTS = function(e) {
            var t, n = "ACGTUNSWMKRYVDHBacgtunswmkryvdhb -",
              r = {};
            for (t = 0; t < n.length; ++t) r[n.charAt(t)] = "TGCAANSWKMYRBHDVtgcaanswkmyrbhdv -".charAt(t);
            return r
          }()), "" === (t = t || "")) return "";
        for (var n = t.split(""), r = 0; r < n.length; ++r) {
          if (!e.COMPLEMENTS[n[r]]) throw "Unknown base: " + n[r];
          n[r] = e.COMPLEMENTS[n[r]]
        }
        for (var i = "", o = n.length - 1; o >= 0; o--) i += n[o];
        return i
      }, r.TmCalc.prototype.calcTm = function(e, t, n) {
        var r, i, o, a, s, u, l, c = 0,
          d = 0,
          p = 0,
          f = 0,
          h = 0,
          g = 0,
          m = 0,
          v = 0,
          $ = 0,
          b = 0,
          y = this.primerc,
          w = 0;
        if (!e && !this.wseq) return {};
        if (t && (!n || "" === n)) throw {
          name: "Missing pairmap",
          message: "All 3 parameters must be present (seq, cseq, pairmap) to calculate a Tm with possible mismatches."
        };
        if (t) this.setSeq(e, t);
        else {
          for (;
            "-" === e[0];) e = e.slice(1);
          for (;
            "-" === e[e.length - 1];) e = e.slice(0, e.length - 1);
          this.setSeq(e)
        }
        if (o = this.wseq.split(""), a = this.cseq.split("").reverse(), n && (this.pairmap = n.toLowerCase()), this.pairmap.length !== this.wseq.length)
          throw {
            name: "Invalid pairmap or alignment",
            message: "Pairmap length " + this.pairmap.length + " must be same as aligned seq length " + this.wseq.length
          };
        if (this.cseq.length !== this.wseq.length) throw {
          name: "Invalid pairmap or alignment",
          message: "aligned seqs must have same length (" + this.wseq.length + ")."
        };
        switch (s = this.pairmap.split(""), this.method) {
          case 1:
            for (f += this.sym ? -12.4 : -10.8, r = 0; r < this.wseq.length - 1; ++r) i = this.wseq.slice(r, r + 2), d += this.dSBr[i], p += this.dHBr[i];
            c = (p *= 1e3) / (f + d + this.R * Math.log(y)) - 273.15 + this.sc_sch;
            break;
          case 3:
          case 4:
            for (g = this.sym ? -1.4 : 0, f = -5.7, h = 200, u = 0, l = this.wseq.length - 1; o[u] + o[u + 1] === "--" || a[u] + a[u + 1] === "--" ||
              "." === s[u];) ++u;
            for ("-" !== o[u] && "-" !== a[u] || (i = o[u] + o[u + 1] + "/" + a[u] + a[u + 1], this.nnde.hasOwnProperty(i) && (m += this.nnde[i].ds, v +=
                this.nnde[i].dh), ++u); o[l - 1] + o[l] === "--" || a[l - 1] + a[l] === "--" || "." === s[l];) --l;
            for ("-" !== o[l] && "-" !== a[l] || (i = o[l - 1] + o[l] + "/" + a[l - 1] + a[l], this.nnde.hasOwnProperty(i) && (m += this.nnde[i].ds, v +=
                this.nnde[i].dh), --l), i = a[u + 1] + a[u] + "/" + o[u + 1] + o[u], this.nntmm.hasOwnProperty(i) && ($ += this.nntmm[i].ds, b += this
                .nntmm[i].dh, ++u), i = o[l - 1] + o[l] + "/" + a[l - 1] + a[l], this.nntmm.hasOwnProperty(i) && ($ += this.nntmm[i].ds, b += this.nntmm[
                i].dh, --l), "a" !== o[u] && "t" !== o[u] || (f += 6.9, h += 2200), "a" !== o[l] && "t" !== o[l] || (f += 6.9, h += 2200), r = u; r <
              l; r++) "|" !== s[r] && "x" !== s[r] || "|" !== s[r + 1] && "x" !== s[r + 1] || (i = o[r] + o[r + 1] + "/" + a[r] + a[r + 1], this.nndhds
              .hasOwnProperty(i) && (d += this.nndhds[i].ds, p += this.nndhds[i].dh));
            var x, C, S;
            w = (c = (1e3 * (p += v + b) + h) / (f + g + (d += m + $) + this.R * Math.log(y))) - 273.15;
            for (var A, T, k, D = /[^_]_+[^_]/g, M = /[^=]=+[^=]/g, E = /[^\^]\^+[^\^]/g, O = 0; A = D.exec(this.pairmap);) T = A.index, x = 2 * ((k = A[
              0].length) - 2), O += this.nntmm[o[T] + o[T + 1] + "/" + a[T] + a[T + 1]].ds, O += this.nntmm[o[T + k - 2] + o[T + k - 1] + "/" + a[T +
              k - 2] + a[T + k - 1]].ds, O += this.nnloop["" + x].ds, ("a" === o[T] && "t" === a[T] || "t" === o[T] && "a" === a[T]) && (O += .5), (
              "a" === o[T + k - 1] && "t" === a[T + k - 1] || "t" === o[T + k - 1] && "a" === a[T + k - 1]) && (O += .5);
            for (; A = M.exec(this.pairmap);) {
              for (T = A.index, x = C = (k = A[0].length) - 2, r = 0; r < k; ++r) "-" === o[T + 1] && (x -= 1), "-" === a[T + 1] && (C -= 1);
              O += this.nntmm[o[T] + o[T + 1] + "/" + a[T] + a[T + 1]].ds, O += this.nntmm[o[T + k - 2] + o[T + k - 1] + "/" + a[T + k - 2] + a[T + k -
                1]].ds, O += .3 * Math.abs(x - C) * 1e3 / 310.15, O += this.nnloop["" + (x + C)].ds, ("a" === o[T] && "t" === a[T] || "t" === o[T] &&
                "a" === a[T]) && (O += .5), ("a" === o[T + k - 1] && "t" === a[T + k - 1] || "t" === o[T + k - 1] && "a" === a[T + k - 1]) && (O += .5)
            }
            for (; A = E.exec(this.pairmap);) T = A.index, 1 == (S = (k = A[0].length) - 2) ? (O += this.nnbulge["" + S].ds, d += this.nndhds[o[T] + o[T +
              2] + "/" + a[T] + a[T + 2]].ds, p += this.nndhds[o[T] + o[T + 2] + "/" + a[T] + a[T + 2]].dh) : S > 1 && (O += this.nnbulge["" + S].ds), (
              "a" === o[T] && "t" === a[T] || "t" === o[T] && "a" === a[T]) && (O += .5), ("a" === o[T + k - 1] && "t" === a[T + k - 1] || "t" === o[T +
              k - 1] && "a" === a[T + k - 1]) && (O += .5);
            w = (c = ((p *= 1e3) + h) / (f + g + (d += O) + this.R * Math.log(y))) - 273.15, 3 === this.method ? c = 1 / (1 / c + this.sc_sl / (p + h)) :
              4 === this.method && (c = 1 / (1 / c + this.sc_ow)), c -= 273.15;
            break;
          case 7:
            c = 81.5 + 16.6 * Math.log(this.saltc) / Math.LN10 + .41 * this.fgc - 675 / this.wseq.length;
            break;
          case 8:
            c = 4 * this.gc_cnt + 2 * this.at_cnt
        }
        return c -= .6 * this.dmso, {
          method: this.method,
          wseq: this.wseq,
          cseq: this.cseq,
          rawtm: w,
          tm: c,
          dh: p + h,
          ds: d + f + g,
          salt: this.saltc,
          ct: this.primerc,
          dmso: this.dmso,
          fgc: this.fgc,
          len: this.wseq.length
        }
      }, r.TmCalc.prototype.calcTm_nomm = function(e) {
        if (!(e && "" !== e || this.wseq)) return {};
        e && (this.setSeq(e), this.saltCorrect());
        var t, n, r = 0,
          i = 0,
          o = 0,
          a = 0,
          s = 0,
          u = 0,
          l = 1e-6 * this.ct,
          c = 0;
        switch (this.method) {
          case 1:
            for (a += this.sym ? -12.4 : -10.8, t = 0; t < this.wseq.length - 1; ++t) n = this.wseq.slice(t, t + 2), i += this.dSBr[n], o += this.dHBr[n];
            r = (c = (o *= 1e3) / (a + i + this.R * Math.log(l)) - 273.15) + this.sc_sch;
            break;
          case 3:
          case 4:
          case 5:
            for (u = this.sym ? -1.4 : 0, "a" !== this.wseq.charAt(0) && "t" !== this.wseq.charAt(0) || (a += 4.1, s += 2300), "a" !== this.wseq.charAt(
                this.wseq.length - 1) && "t" !== this.wseq.charAt(this.wseq.length - 1) || (a += 4.1, s += 2300), "g" !== this.wseq.charAt(0) && "c" !==
              this.wseq.charAt(0) || (a += -2.8, s += 100), "g" !== this.wseq.charAt(this.wseq.length - 1) && "c" !== this.wseq.charAt(this.wseq.length -
                1) || (a += -2.8, s += 100), t = 0; t < this.wseq.length - 1; t++) n = this.wseq.slice(t, t + 2), i += this.dSSa[n], o += this.dHSa[n];
            c = (r = ((o *= 1e3) + s) / (a + u + i + this.R * Math.log(l))) - 273.15, 3 === this.method ? r = 1 / (1 / r + this.sc_sl / (o + s)) : 4 ===
              this.method ? r = 1 / (1 / r + this.sc_ow) : 5 === this.method && (r += this.sc_sch), r -= 273.15;
            break;
          case 7:
            r = 81.5 + 16.6 * Math.log(this.saltc) / Math.LN10 + .41 * this.fgc - 675 / this.wseq.length
        }
        return r -= .6 * this.dmso, {
          method: this.method,
          wseq: this.wseq,
          rawtm: c,
          tm: r,
          dh: o + s,
          ds: i + a + u,
          salt: this.saltc,
          ct: this.primerc,
          dmso: this.dmso,
          fgc: this.fgc,
          len: this.wseq.length
        }
      }, void 0 !== n ? (void 0 !== t && t.exports && (n = t.exports = i), n.NEB_tmelt = i) : (void 0).NEB_tmelt = i, "function" == typeof define &&
      define.amd && define(i, ["NEB_tmelt"], (function() {
        return i
      }))
  }, {}],
  15: [function(e, t, n) {
    "use strict";
    t.exports = function() {
      var e = e || {};
      return e.trim = function(e) {
        for (var t = (e = e.replace(/^\s+/, "")).length - 1; t >= 0; t--)
          if (/\S/.test(e.charAt(t))) {
            e = e.substring(0, t + 1);
            break
          } return e
      }, e.revcomp = function(e) {
        var t, n, r = e.split("").reverse();
        for (t = 0; t < r.length; ++t) {
          if (-1 === (n = "ACGTUKMSWRYBDHVNacgtukmswrybdhvn".indexOf(r[t]))) throw {
            name: "InvalidCharacter",
            message: "A non IUPAC base was encountered"
          };
          r[t] = "TGCAAMKWSYRVHDBNtgcaamkwsyrvhdbn".charAt(n)
        }
        return r.join("")
      }, e.isValidSeq = function(e, t) {
        var n = " acgtu";
        n += n.toUpperCase();
        var r, i, o = " acgtuwsrymkbdhvn";
        for (o += o.toUpperCase(), r = !1 === t ? n : o, i = 0; i < e.length; ++i)
          if (-1 === r.indexOf(e.charAt(i))) return !1;
        return !0
      }, e.nonCanonicalCount = function(e) {
        var t = e.replace(/\s+/, "").split(""),
          n = "wsrymkbdhvn".split("");
        return t.reduce((function(e, t) {
          return -1 != n.indexOf(t) ? e + 1 : e
        }), 0)
      }, e.isNumeric = function(e) {
        return !isNaN(e / 1)
      }, e.dmsg = function(e) {
        "undefined" != typeof console && "function" == typeof console.log && console.log(e)
      }, e
    }
  }, {}],
  16: [function(e, t, n) {
    "use strict";
    t.exports = function(e) {
      return {
        getData: function(t) {
          return t || (t = "data/tmcalculatordata.json"), e.get(t).then((function(e) {
            return e.data
          }))
        }
      }
    }
  }, {}],
  17: [function(e, t, n) {
    "use strict";
    var r = e("../libs/NEB_bioseq.js");
    t.exports = function() {
      return r
    }
  }, {
    "../libs/NEB_bioseq.js": 13
  }],
  18: [function(e, t, n) {
    "use strict";
    var r = e("angular").module("tmcApp");
    r.factory("modalService", ["$http", "$uibModal", e("./modalService2")]), r.factory("nebutil", [e("../libs/neb_util")]), r.factory("tmData", ["$http",
      e("./TmData")
    ]), r.factory("neb_bioseq", [e("./bioseq")]), r.factory("tmcalculatorData", ["$http", "nebutil", e("./tmcalculator")]), r.factory("tmcalc", [e(
      "./tmcalc")])
  }, {
    "../libs/neb_util": 15,
    "./TmData": 16,
    "./bioseq": 17,
    "./modalService2": 19,
    "./tmcalc": 20,
    "./tmcalculator": 21,
    angular: 32
  }],
  19: [function(e, t, n) {
    "use strict";
    t.exports = function(e, t) {
      var n = {
          backdrop: !0,
          keyboard: !0,
          modalFade: !0
        },
        r = {},
        i = {
          showModal: function(e, t) {
            return e || (e = {}), e.backdrop = !0, this.show(e, t)
          },
          show: function(e, i) {
            var o = {},
              a = {};
            return angular.extend(o, n, e), angular.extend(a, r, i), o.controller || (o.controller = ["$scope", "$uibModalInstance", function(e, t) {
              e.modalOptions = a, e.modalOptions.close = function(e) {
                t.close("done")
              }, e.modalOptions.ok = function(e) {
                t.close("done")
              }
            }]), t.open(o).result
          }
        };
      return i
    }
  }, {}],
  20: [function(e, t, n) {
    "use strict";
    var r = e("../libs/NEB_tmelt.js");
    t.exports = function() {
      return new r.TmCalc
    }
  }, {
    "../libs/NEB_tmelt.js": 14
  }],
  21: [function(e, t, n) {
    /*
     * tmcalculator
     * user/repo
     * @preserve
     * Copyright (c) 2015-2018 New England Biolabs
     * Author: Sanjay Kumar
     * Licensed under the AGPLv3 license.
     */
    "use strict";
    t.exports = function(e, t) {
      var n = t;
      return {
        saveddata: {},
        tmcdata: null,
        setData: function(e) {
          this.tmcdata = e
        },
        getData: function() {
          return this.tmcdata
        },
        getGroups: function() {
          return this.tmcdata.groups
        },
        getGroupKeys: function() {
          var e, t = [];
          for (e in this.tmcdata.groups) this.tmcdata.groups.hasOwnProperty(e) && t.push(e);
          return t
        },
        getProducts: function() {
          return this.tmcdata.prods
        },
        getBuffers: function() {
          return this.tmcdata.buffs
        },
        getBufferSaltForProduct: function(e) {
          return this.tmcdata.buffs[this.tmcdata.prods[e].buffer]
        },
        getBufferIdForProduct: function(e) {
          return this.tmcdata.prods[e].buffer
        },
        getProductsForGroup: function(e) {
          for (var t, n = this.tmcdata.groups[e], r = [], i = 0; i < n.length; ++i)(t = this.tmcdata.prods[n[i]]).id = n[i], r.push(t);
          return r
        },
        getProductKeysForGroup: function(e) {
          var t, n = [],
            r = this.tmcdata.groups[e];
          for (t in r) n.push(this.tmcdata.prods[r[t]].name);
          return n
        },
        getPropsForProduct: function(e) {
          return this.tmcdata.prods[e]
        },
        saveUserPrefs: function(e) {
          this.saveddata = e
        },
        restoreUserPrefs: function() {
          return this.savedata
        },
        validateInput: function(e, t, r, i, o) {
          var a = [],
            s = !1,
            u = [],
            l = e.length,
            c = t.length,
            d = !0,
            p = !0,
            f = !0,
            h = "",
            g = "";
          return (r <= 0 || isNaN(r / 1)) && (d = !1, s = !0, a.push("Invalid primer concentration. ")), r < .5 && -1 !== o.indexOf("Q5") ? u.push(
              "The recommended primer concentration for Q5/Q5U reactions is 500 nM.") : r < .5 && -1 !== o.indexOf("Phusion") ? u.push(
              "The recommended primer concentration for Phusion reactions is 500 nM.") : r < .4 && -1 !== o.indexOf("LongAmp") ? u.push(
              "The recommended primer concentration for LongAmp reactions is 400 nM.") : r < .2 && u.push(
              "The recommended primer concentration is 200 nM."), p = n.isValidSeq(e, !0), f = n.isValidSeq(t, !0), p = p && n.nonCanonicalCount(e) <
            4, f = f && n.nonCanonicalCount(t) < 4, p || f ? (p || (a.push("Primer 1 has invalid bases or more than 3 ambiguous bases. "), s = !0, h =
              "invalidseq"), f || (a.push("Primer 2 has invalid bases or more than 3 ambiguous bases. "), s = !0, g = "invalidseq")) : (a.push(
              "Both primers have invalid bases or more than 3 ambiguous bases."), s = !0, h = "invalidseq", g = "invalidseq"), 0 === l && (a.push(
              "Primer 1 missing. "), s = !0, h = "invalidseq"), 0 === c && (a.push("Primer 2 missing. "), s = !0, g = "invalidseq"), (l < 8 || c <
            8) && ((l > 0 || c > 0) && (a.push("Both primers need to be longer than 7 nt"), s = !0), l < 8 && (h = "invalidseq"), c < 8 && (g =
              "invalidseq")), {
              hasWarnings: !1,
              warnings: u,
              hasCritWarnings: s,
              critwarnings: a,
              ctisValid: d,
              p1isValid: p,
              p2isValid: f,
              p1status: h,
              p2status: g
            }
        },
        validateTm: function(e, t, n, r, i, o) {
          var a = [],
            s = [];
          return (e - t) * (e - t) > 25 && s.push("Tm difference is greater than the recommended limit of 5 °C. "), "" !== n && n < 60 && -1 !== i
            .indexOf("Q5U") ? s.push(" The minimum recommended annealing temperature for Q5U is 60 °C.") : "" !== n && n < 55 && -1 !== i.indexOf(
              "Q5") ? s.push(" The minimum recommended annealing temperature for Q5 is 55 °C.") : "" !== n && n < 45 && s.push(
              " Annealing temperature is lower than the recommended minimum of 45 °C."), "" !== n && n >= 65 && (n < e || n < t) && (72 !== n ||
              "Phusion" !== i && "Vent" !== i && "Deep Vent" !== i && "Phusion Flex" !== i && -1 === i.indexOf("Q5") ? 65 !== n || "LongAmp Taq" !==
              i && "LongAmp Hot Start Taq" !== i ? 68 === n && a.push(
                "Annealing temperature for experiments with this enzyme should typically not exceed 68°C.") : a.push(
                "Annealing temperature for experiments with this enzyme should typically not exceed 65°C.") : a.push(
                "Annealing temperature for experiments with this enzyme should typically not exceed 72°C.")), {
              hasWarnings: a.length > 0,
              warnings: a,
              hasCritWarnings: s.length > 0,
              critwarnings: s
            }
        },
        validateBuffer: function(e, t, n) {
          var r = [],
            i = [];
          return "phusion_gc" === n || "onetaq_gc" === n || "phusionflex_gc" === n ? r.push(
            "DMSO can improve PCR amplification from GC-rich templates,  but it is also known to reduce the annealing temperature of primers in a PCR reaction. Therefore,  it is recommended that for every 1% of additional DMSO added, the calculated annealing temperature should be reduced by 0.6°C [Chester and Marshak,  1993. Analytical Biochemistry 209,  284-290]."
            ) : "q5" === n ? r.push(
            "Use of the Q5 High GC Enhancer often lowers the range of temperatures at which specific amplification can be observed, however the rule used to determine Q5 annealing temperatures (Ta = Tm_lower+1°C) typically yields values that will support specific amplification with or without the enhancer."
            ) : "q5u" === n && r.push("For bisulfite-converted or deaminated DNA substrates, a minimum Ta of 60°C should be used."), {
            hasWarnings: r.length > 0,
            warnings: r,
            hasCritWarnings: i.length > 0,
            critwarnings: i
          }
        },
        getAnnealTemp: function(e, t, n, r, i, o) {
          console.log("getAnnealTemp"), console.log(i, o), o = o.toLowerCase();
          var a, s, u, l = e.replace(/\s/g, "").length,
            c = n.replace(/\s/g, "").length;
          switch (u = l < c ? l : c, a = s = t < r ? t : r, i) {
            case "Taq DNA Polymerase":
            case "Hemo KlenTaq":
            case "OneTaq":
            case "OneTaq Hot Start":
            case "EpiMark Hot Start":
            case "Hot Start Taq":
              u > 7 && (a = s - 5), a > 68 && (a = 68);
              break;
            case "LongAmp Taq":
            case "LongAmp Hot Start Taq":
              u > 7 && (a = s - 5), a > 65 && (a = 65);
              break;
            case "Vent":
            case "Deep Vent":
              u > 20 && (a = s - 2), a > 72 && (a = 72);
              break;
            case "Phusion":
            case "Phusion Hot Start Flex":
              console.log("minlen", u), (a = .93 * s + 7.5) > 72 && (a = 72);
              break;
            case "Q5U Hot Start":
              u > 7 && (a = s + 2), a > 72 && (a = 72);
              break;
            case "Q5":
            case "Q5 Hot Start":
            case "Q5 Blood Direct":
              u > 7 && (a = s + 1), a > 72 && (a = 72);
              break;
            case "Master Mix":
              0 === o.indexOf("q5u") ? (u > 7 && (a = s + 2), a > 72 && (a = 72)) : 0 === o.indexOf("q5") ? (u > 7 && (a = s + 1), a > 72 && (a =
                72)) : 0 === o.indexOf("phusion") ? (a = .93 * s + 7.5) > 72 && (a = 72) : 0 === o.indexOf("vent") || 0 === o.indexOf("deep") ? (u >
                  20 && (a = s - 2), a > 72 && (a = 72)) : (u > 7 && (a = s - 5), a > 68 && (a = 68));
              break;
            default:
              console.log("using default case as group " + i + " not matched"), a = s
          }
          return Math.round(10 * a) / 10
        }
      }
    }
  }, {}],
  22: [function(e, t, n) {
    "use strict";
    var r = r || function(e) {
      if ("undefined" == typeof navigator || !/MSIE [1-9]\./.test(navigator.userAgent)) {
        var t = e.document,
          n = function() {
            return e.URL || e.webkitURL || e
          },
          r = t.createElementNS("http://www.w3.org/1999/xhtml", "a"),
          i = "download" in r,
          o = e.webkitRequestFileSystem,
          a = e.requestFileSystem || o || e.mozRequestFileSystem,
          s = function(t) {
            (e.setImmediate || e.setTimeout)((function() {
              throw t
            }), 0)
          },
          u = "application/octet-stream",
          l = 0,
          c = function(t) {
            var r = function() {
              "string" == typeof t ? n().revokeObjectURL(t) : t.remove()
            };
            e.chrome ? r() : setTimeout(r, 500)
          },
          d = function(e, t, n) {
            for (var r = (t = [].concat(t)).length; r--;) {
              var i = e["on" + t[r]];
              if ("function" == typeof i) try {
                i.call(e, n || e)
              } catch (e) {
                s(e)
              }
            }
          },
          p = function(e) {
            return /^\s*(?:text\/\S*|application\/xml|\S*\/\S*\+xml)\s*;.*charset\s*=\s*utf-8/i.test(e.type) ? new Blob(["\ufeff", e], {
              type: e.type
            }) : e
          },
          f = function(t, s, f) {
            f || (t = p(t));
            var h, g, m, v = this,
              $ = t.type,
              b = !1,
              y = function() {
                d(v, "writestart progress write writeend".split(" "))
              },
              w = function() {
                (!b && h || (h = n().createObjectURL(t)), g) ? g.location.href = h: null == e.open(h, "_blank") && "undefined" != typeof safari && (e
                  .location.href = h);
                v.readyState = v.DONE, y(), c(h)
              },
              x = function(e) {
                return function() {
                  if (v.readyState !== v.DONE) return e.apply(this, arguments)
                }
              },
              C = {
                create: !0,
                exclusive: !1
              };
            if (v.readyState = v.INIT, s || (s = "download"), i) return h = n().createObjectURL(t), r.href = h, r.download = s, void setTimeout((
              function() {
                var e, t;
                e = r, t = new MouseEvent("click"), e.dispatchEvent(t), y(), c(h), v.readyState = v.DONE
              }));
            e.chrome && $ && $ !== u && (m = t.slice || t.webkitSlice, t = m.call(t, 0, t.size, u), b = !0), o && "download" !== s && (s +=
              ".download"), ($ === u || o) && (g = e), a ? (l += t.size, a(e.TEMPORARY, l, x((function(e) {
                e.root.getDirectory("saved", C, x((function(e) {
                  var n = function() {
                    e.getFile(s, C, x((function(e) {
                      e.createWriter(x((function(n) {
                        n.onwriteend = function(t) {
                          g.location.href = e.toURL(), v.readyState = v.DONE, d(v, "writeend", t), c(e)
                        }, n.onerror = function() {
                          var e = n.error;
                          e.code !== e.ABORT_ERR && w()
                        }, "writestart progress write abort".split(" ").forEach((function(e) {
                          n["on" + e] = v["on" + e]
                        })), n.write(t), v.abort = function() {
                          n.abort(), v.readyState = v.DONE
                        }, v.readyState = v.WRITING
                      })), w)
                    })), w)
                  };
                  e.getFile(s, {
                    create: !1
                  }, x((function(e) {
                    e.remove(), n()
                  })), x((function(e) {
                    e.code === e.NOT_FOUND_ERR ? n() : w()
                  })))
                })), w)
              })), w)) : w()
          },
          h = f.prototype;
        return "undefined" != typeof navigator && navigator.msSaveOrOpenBlob ? function(e, t, n) {
          return n || (e = p(e)), navigator.msSaveOrOpenBlob(e, t || "download")
        } : (h.abort = function() {
            var e = this;
            e.readyState = e.DONE, d(e, "abort")
          }, h.readyState = h.INIT = 0, h.WRITING = 1, h.DONE = 2, h.error = h.onwritestart = h.onprogress = h.onwrite = h.onabort = h.onerror = h
          .onwriteend = null,
          function(e, t, n) {
            return new f(e, t, n)
          })
      }
    }("undefined" != typeof self && self || "undefined" != typeof window && window || (void 0).content);
    void 0 !== t && t.exports ? t.exports.saveAs = r : "undefined" != typeof define && null !== define && null != define.amd && define([], (function() {
      return r
    }))
  }, {}],
  23: [function(e, t, n) {
    /**
     * @license AngularJS v1.8.2
     * (c) 2010-2020 Google LLC. http://angularjs.org
     * License: MIT
     */
    ! function(e, t) {
      "use strict";
      var n, r, i, o, a = "-add",
        s = "-remove",
        u = "ng-animate",
        l = "$$ngAnimateChildren";
      void 0 === e.ontransitionend && void 0 !== e.onwebkittransitionend ? ("-webkit-", n = "WebkitTransition", r = "webkitTransitionEnd transitionend") :
        (n = "transition", r = "transitionend"), void 0 === e.onanimationend && void 0 !== e.onwebkitanimationend ? ("-webkit-", i = "WebkitAnimation",
          o = "webkitAnimationEnd animationend") : (i = "animation", o = "animationend");
      var c = "Duration",
        d = "Property",
        p = "Delay",
        f = "TimingFunction",
        h = i + p,
        g = i + c,
        m = n + p,
        v = n + c,
        $ = t.$$minErr("ng");

      function b(e, t, n) {
        if (!e) throw $("areq", "Argument '{0}' is {1}", t || "?", n || "required");
        return e
      }

      function y(e, t) {
        return e || t ? e ? t ? (Y(e) && (e = e.join(" ")), Y(t) && (t = t.join(" ")), e + " " + t) : e : t : ""
      }

      function w(e, t, n) {
        var r = "";
        return e = Y(e) ? e : e && Z(e) && e.length ? e.split(/\s+/) : [], z(e, (function(e, i) {
          e && e.length > 0 && (r += i > 0 ? " " : "", r += n ? t + e : e + t)
        })), r
      }

      function x(e) {
        if (e instanceof te) switch (e.length) {
          case 0:
            return e;
          case 1:
            if (1 === e[0].nodeType) return e;
            break;
          default:
            return te(C(e))
        }
        if (1 === e.nodeType) return te(e)
      }

      function C(e) {
        if (!e[0]) return e;
        for (var t = 0; t < e.length; t++) {
          var n = e[t];
          if (1 === n.nodeType) return n
        }
      }

      function S(e) {
        return function(t, n) {
          n.addClass && (! function(e, t, n) {
            z(t, (function(t) {
              e.addClass(t, n)
            }))
          }(e, t, n.addClass), n.addClass = null), n.removeClass && (! function(e, t, n) {
            z(t, (function(t) {
              e.removeClass(t, n)
            }))
          }(e, t, n.removeClass), n.removeClass = null)
        }
      }

      function A(e) {
        if (!(e = e || {}).$$prepared) {
          var t = e.domOperation || ne;
          e.domOperation = function() {
            e.$$domOperationFired = !0, t(), t = ne
          }, e.$$prepared = !0
        }
        return e
      }

      function T(e, t) {
        k(e, t), D(e, t)
      }

      function k(e, t) {
        t.from && (e.css(t.from), t.from = null)
      }

      function D(e, t) {
        t.to && (e.css(t.to), t.to = null)
      }

      function M(e, t, n) {
        var r = t.options || {},
          i = n.options || {},
          o = (r.addClass || "") + " " + (i.addClass || ""),
          u = (r.removeClass || "") + " " + (i.removeClass || ""),
          l = function(e, t, n) {
            var r = 1,
              i = -1,
              o = {};
            e = l(e), t = l(t), z(t, (function(e, t) {
              o[t] = r
            })), n = l(n), z(n, (function(e, t) {
              o[t] = o[t] === r ? null : i
            }));
            var u = {
              addClass: "",
              removeClass: ""
            };

            function l(e) {
              Z(e) && (e = e.split(" "));
              var t = {};
              return z(e, (function(e) {
                e.length && (t[e] = !0)
              })), t
            }
            return z(o, (function(t, n) {
              var o, l;
              t === r ? (o = "addClass", l = !e[n] || e[n + s]) : t === i && (o = "removeClass", l = e[n] || e[n + a]), l && (u[o].length && (u[
                o] += " "), u[o] += n)
            })), u
          }(e.attr("class"), o, u);
        i.preparationClasses && (r.preparationClasses = P(i.preparationClasses, r.preparationClasses), delete i.preparationClasses);
        var c = r.domOperation !== ne ? r.domOperation : null;
        return W(r, i), c && (r.domOperation = c), l.addClass ? r.addClass = l.addClass : r.addClass = null, l.removeClass ? r.removeClass = l
          .removeClass : r.removeClass = null, t.addClass = r.addClass, t.removeClass = r.removeClass, r
      }

      function E(e) {
        return e instanceof te ? e[0] : e
      }

      function O(e, t) {
        var n = t ? "paused" : "",
          r = i + "PlayState";
        return N(e, [r, n]), [r, n]
      }

      function N(e, t) {
        var n = t[0],
          r = t[1];
        e.style[n] = r
      }

      function P(e, t) {
        return e ? t ? e + " " + t : e : t
      }
      var q = function(e, t) {
          var n = t ? "-" + t + "s" : "";
          return N(e, [m, n]), [m, n]
        },
        I = ["$interpolate", function(e) {
          return {
            link: function(t, n, r) {
              var i = r.ngAnimateChildren;

              function o(e) {
                e = "on" === e || "true" === e, n.data(l, e)
              }
              Z(i) && 0 === i.length ? n.data(l, !0) : (o(e(i)(t)), r.$observe("ngAnimateChildren", o))
            }
          }
        }],
        L = "$$animateCss",
        R = 1e3,
        j = {
          transitionDuration: v,
          transitionDelay: m,
          transitionProperty: n + d,
          animationDuration: g,
          animationDelay: h,
          animationIterationCount: i + "IterationCount"
        },
        _ = {
          transitionDuration: v,
          transitionDelay: m,
          animationDuration: g,
          animationDelay: h
        };

      function V(e, t) {
        return [t ? h : m, e + "s"]
      }

      function U(e, t, n) {
        var r = Object.create(null),
          i = e.getComputedStyle(t) || {};
        return z(n, (function(e, t) {
          var n, o, a = i[e];
          if (a) {
            var s = a.charAt(0);
            ("-" === s || "+" === s || s >= 0) && (n = 0, o = a.split(/\s*,\s*/), z(o, (function(e) {
              "s" === e.charAt(e.length - 1) && (e = e.substring(0, e.length - 1)), e = parseFloat(e) || 0, n = n ? Math.max(e, n) : e
            })), a = n), 0 === a && (a = null), r[t] = a
          }
        })), r
      }

      function H(e) {
        return 0 === e || null != e
      }

      function F(e, t) {
        var r = n,
          i = e + "s";
        return t ? r += c : i += " linear all", [r, i]
      }

      function B(e, t, n) {
        z(n, (function(n) {
          e[n] = K(e[n]) ? e[n] : t.style.getPropertyValue(n)
        }))
      }
      var G, W, z, Y, K, Q, X, J, Z, ee, te, ne, re = ["$animateProvider", function(e) {
          this.$get = ["$window", "$$jqLite", "$$AnimateRunner", "$timeout", "$$animateCache", "$$forceReflow", "$sniffer", "$$rAFScheduler",
            "$$animateQueue",
            function(e, t, u, l, c, p, h, m, v) {
              var $ = S(t);
              var b = [];

              function y(e) {
                b.push(e), m.waitUntilQuiet((function() {
                  c.flush();
                  for (var e = p(), t = 0; t < b.length; t++) b[t](e);
                  b.length = 0
                }))
              }

              function x(t, n, r, i) {
                var o = function(t, n, r, i, o) {
                    var a = c.get(r);
                    a || "infinite" === (a = U(e, t, o)).animationIterationCount && (a.animationIterationCount = 1);
                    var s = i || a.transitionDuration > 0 || a.animationDuration > 0;
                    return c.put(r, a, s), a
                  }(t, 0, r, i, j),
                  a = o.animationDelay,
                  s = o.transitionDelay;
                return o.maxDelay = a && s ? Math.max(a, s) : a || s, o.maxDuration = Math.max(o.animationDuration * o.animationIterationCount, o
                  .transitionDuration), o
              }
              return function(p, m) {
                var b = m || {};
                b.$$prepared || (b = A(G(b)));
                var C = {},
                  S = E(p);
                if (!S || !S.parentNode || !v.enabled()) return De();
                var M, P, I, j, W, K, Q, X, J, Z, ee = [],
                  te = (p.attr("class"), function(e) {
                    var t = {};
                    return e && (e.to || e.from) && (t.to = e.to, t.from = e.from), t
                  }(b)),
                  re = [];
                if (0 === b.duration || !h.animations && !h.transitions) return De();
                var ie = b.event && Y(b.event) ? b.event.join(" ") : b.event,
                  oe = ie && b.structural,
                  ae = "",
                  se = "";
                oe ? ae = w(ie, "ng-", !0) : ie && (ae = ie), b.addClass && (se += w(b.addClass, a)), b.removeClass && (se.length && (se += " "),
                  se += w(b.removeClass, s)), b.applyClassesEarly && se.length && $(p, b);
                var ue = [ae, se].join(" ").trim(),
                  le = te.to && Object.keys(te.to).length > 0;
                if (!((b.keyframeStyle || "").length > 0) && !le && !ue) return De();
                var ce, de, pe = c.cacheKey(S, ie, b.addClass, b.removeClass);
                if (c.containsCachedAnimationWithoutDuration(pe)) return ue = null, De();
                if (b.stagger > 0) {
                  var fe = parseFloat(b.stagger);
                  ce = {
                    transitionDelay: fe,
                    animationDelay: fe,
                    transitionDuration: 0,
                    animationDuration: 0
                  }
                } else ce = function(n, r, i, o) {
                  var a, s = "stagger-" + i;
                  if (c.count(i) > 0 && !(a = c.get(s))) {
                    var u = w(r, "-stagger");
                    t.addClass(n, u), (a = U(e, n, o)).animationDuration = Math.max(a.animationDuration, 0), a.transitionDuration = Math.max(a
                      .transitionDuration, 0), t.removeClass(n, u), c.put(s, a, !0)
                  }
                  return a || {}
                }(S, ue, pe, _);
                if (b.$$skipPreparationClasses || t.addClass(p, ue), b.transitionStyle) {
                  var he = [n, b.transitionStyle];
                  N(S, he), ee.push(he)
                }
                if (b.duration >= 0) {
                  de = S.style[n].length > 0;
                  var ge = F(b.duration, de);
                  N(S, ge), ee.push(ge)
                }
                if (b.keyframeStyle) {
                  var me = [i, b.keyframeStyle];
                  N(S, me), ee.push(me)
                }
                var ve = ce ? b.staggerIndex >= 0 ? b.staggerIndex : c.count(pe) : 0,
                  $e = 0 === ve;
                $e && !b.skipBlocking && q(S, 9999);
                var be = x(S, 0, pe, !oe),
                  ye = be.maxDelay;
                K = Math.max(ye, 0), X = be.maxDuration;
                var we = {};
                if (we.hasTransitions = be.transitionDuration > 0, we.hasAnimations = be.animationDuration > 0, we.hasTransitionAll = we
                  .hasTransitions && "all" === be.transitionProperty, we.applyTransitionDuration = le && (we.hasTransitions && !we
                    .hasTransitionAll || we.hasAnimations && !we.hasTransitions), we.applyAnimationDuration = b.duration && we.hasAnimations, we
                  .applyTransitionDelay = H(b.delay) && (we.applyTransitionDuration || we.hasTransitions), we.applyAnimationDelay = H(b.delay) && we
                  .hasAnimations, we.recalculateTimingStyles = se.length > 0, (we.applyTransitionDuration || we.applyAnimationDuration) && (X = b
                    .duration ? parseFloat(b.duration) : X, we.applyTransitionDuration && (we.hasTransitions = !0, be.transitionDuration = X, de = S
                      .style[n + d].length > 0, ee.push(F(X, de))), we.applyAnimationDuration && (we.hasAnimations = !0, be.animationDuration = X,
                      ee.push([g, X + "s"]))), 0 === X && !we.recalculateTimingStyles) return De();
                var xe, Ce = w(ue, "-active");
                null != b.delay && ("boolean" != typeof b.delay && (xe = parseFloat(b.delay), K = Math.max(xe, 0)), we.applyTransitionDelay && ee
                  .push(V(xe)), we.applyAnimationDelay && ee.push(V(xe, !0)));
                return null == b.duration && be.transitionDuration > 0 && (we.recalculateTimingStyles = we.recalculateTimingStyles || $e), Q = K *
                  R, J = X * R, b.skipBlocking || (we.blockTransition = be.transitionDuration > 0, we.blockKeyframeAnimation = be
                    .animationDuration > 0 && ce.animationDelay > 0 && 0 === ce.animationDuration), b.from && (b.cleanupStyles && B(C, S, Object
                    .keys(b.from)), k(p, b)), we.blockTransition || we.blockKeyframeAnimation ? ke(X) : b.skipBlocking || q(S, !1), {
                    $$willAnimate: !0,
                    end: Se,
                    start: function() {
                      if (!M) return j = new u(W = {
                        end: Se,
                        cancel: Ae,
                        resume: null,
                        pause: null
                      }), y(Ee), j
                    }
                  };

                function Se() {
                  Te()
                }

                function Ae() {
                  Te(!0)
                }

                function Te(e) {
                  if (!(M || I && P)) {
                    M = !0, P = !1, ue && !b.$$skipPreparationClasses && t.removeClass(p, ue), Ce && t.removeClass(p, Ce), O(S, !1), q(S, !1), z(ee,
                      (function(e) {
                        S.style[e[0]] = ""
                      })), $(p, b), T(p, b), Object.keys(C).length && z(C, (function(e, t) {
                      e ? S.style.setProperty(t, e) : S.style.removeProperty(t)
                    })), b.onDone && b.onDone(), re && re.length && p.off(re.join(" "), Me);
                    var n = p.data(L);
                    n && (l.cancel(n[0].timer), p.removeData(L)), j && j.complete(!e)
                  }
                }

                function ke(e) {
                  we.blockTransition && q(S, e), we.blockKeyframeAnimation && O(S, !!e)
                }

                function De() {
                  return j = new u({
                    end: Se,
                    cancel: Ae
                  }), y(ne), Te(), {
                    $$willAnimate: !1,
                    start: function() {
                      return j
                    },
                    end: Se
                  }
                }

                function Me(e) {
                  e.stopPropagation();
                  var t = e.originalEvent || e;
                  if (t.target === S) {
                    var n = t.$manualTimeStamp || Date.now(),
                      r = parseFloat(t.elapsedTime.toFixed(3));
                    Math.max(n - Z, 0) >= Q && r >= X && (I = !0, Te())
                  }
                }

                function Ee() {
                  if (!M)
                    if (S.parentNode) {
                      var e = function(e) {
                          if (I) P && e && (P = !1, Te());
                          else if (P = !e, be.animationDuration) {
                            var t = O(S, P);
                            P ? ee.push(t) : (r = t, i = (n = ee).indexOf(r), r >= 0 && n.splice(i, 1))
                          }
                          var n, r, i
                        },
                        a = ve > 0 && (be.transitionDuration && 0 === ce.transitionDuration || be.animationDuration && 0 === ce
                        .animationDuration) && Math.max(ce.animationDelay, ce.transitionDelay);
                      a ? l(s, Math.floor(a * ve * R), !1) : s(), W.resume = function() {
                        e(!0)
                      }, W.pause = function() {
                        e(!1)
                      }
                    } else Te();

                  function s() {
                    if (!M) {
                      if (ke(!1), z(ee, (function(e) {
                          var t = e[0],
                            n = e[1];
                          S.style[t] = n
                        })), $(p, b), t.addClass(p, Ce), we.recalculateTimingStyles) {
                        if (S.getAttribute("class") + " " + ue, pe = c.cacheKey(S, ie, b.addClass, b.removeClass), be = x(S, 0, pe, !1), ye = be
                          .maxDelay, K = Math.max(ye, 0), 0 === (X = be.maxDuration)) return void Te();
                        we.hasTransitions = be.transitionDuration > 0, we.hasAnimations = be.animationDuration > 0
                      }
                      if (we.applyAnimationDelay && (ye = "boolean" != typeof b.delay && H(b.delay) ? parseFloat(b.delay) : ye, K = Math.max(ye, 0),
                          be.animationDelay = ye, xe = V(ye, !0), ee.push(xe), S.style[xe[0]] = xe[1]), Q = K * R, J = X * R, b.easing) {
                        var e, a = b.easing;
                        we.hasTransitions && (e = n + f, ee.push([e, a]), S.style[e] = a), we.hasAnimations && (e = i + f, ee.push([e, a]), S.style[
                          e] = a)
                      }
                      be.transitionDuration && re.push(r), be.animationDuration && re.push(o), Z = Date.now();
                      var s = Q + 1.5 * J,
                        d = Z + s,
                        h = p.data(L) || [],
                        g = !0;
                      if (h.length) {
                        var m = h[0];
                        (g = d > m.expectedEndTime) ? l.cancel(m.timer): h.push(Te)
                      }
                      if (g) {
                        var v = l(u, s, !1);
                        h[0] = {
                          timer: v,
                          expectedEndTime: d
                        }, h.push(Te), p.data(L, h)
                      }
                      re.length && p.on(re.join(" "), Me), b.to && (b.cleanupStyles && B(C, S, Object.keys(b.to)), D(p, b))
                    }
                  }

                  function u() {
                    var e = p.data(L);
                    if (e) {
                      for (var t = 1; t < e.length; t++) e[t]();
                      p.removeData(L)
                    }
                  }
                }
              }
            }
          ]
        }],
        ie = ["$$animationProvider", function(e) {
          e.drivers.push("$$animateCssDriver");
          var t = "ng-animate-shim",
            n = "ng-anchor",
            r = "ng-anchor-out";
          this.$get = ["$animateCss", "$rootScope", "$$AnimateRunner", "$rootElement", "$sniffer", "$$jqLite", "$document", function(e, i, o, a, s, u,
            l) {
            if (!s.animations && !s.transitions) return ne;
            var c, d = l[0].body,
              p = E(a),
              f = te((c = p).parentNode && 11 === c.parentNode.nodeType || d.contains(p) ? p : d);
            return function(i) {
              return i.from && i.to ? function(i, a, s, u) {
                var l = m(i),
                  c = m(a),
                  p = [];
                if (z(u, (function(i) {
                    var a = function(i, a, s) {
                      var u = te(E(a).cloneNode(!0)),
                        l = h(b(u));
                      a.addClass(t), s.addClass(t), u.addClass(n), f.append(u);
                      var c, p = $();
                      if (!p && !(c = y())) return w();
                      var m = p || c;
                      return {
                        start: function() {
                          var e, t = m.start();
                          return t.done((function() {
                            if (t = null, !c && (c = y())) return (t = c.start()).done((function() {
                              t = null, w(), e.complete()
                            })), t;
                            w(), e.complete()
                          })), e = new o({
                            end: n,
                            cancel: n
                          });

                          function n() {
                            t && t.end()
                          }
                        }
                      };

                      function v(e) {
                        var t = {},
                          n = E(e).getBoundingClientRect();
                        return z(["width", "height", "top", "left"], (function(e) {
                          var r = n[e];
                          switch (e) {
                            case "top":
                              r += d.scrollTop;
                              break;
                            case "left":
                              r += d.scrollLeft
                          }
                          t[e] = Math.floor(r) + "px"
                        })), t
                      }

                      function $() {
                        var t = e(u, {
                          addClass: r,
                          delay: !0,
                          from: v(a)
                        });
                        return t.$$willAnimate ? t : null
                      }

                      function b(e) {
                        return e.attr("class") || ""
                      }

                      function y() {
                        var t = h(b(s)),
                          n = g(t, l),
                          r = g(l, t),
                          i = e(u, {
                            to: v(s),
                            addClass: "ng-anchor-in " + n,
                            removeClass: "ng-anchor-out " + r,
                            delay: !0
                          });
                        return i.$$willAnimate ? i : null
                      }

                      function w() {
                        u.remove(), a.removeClass(t), s.removeClass(t)
                      }
                    }(0, i.out, i.in);
                    a && p.push(a)
                  })), !l && !c && 0 === p.length) return;
                return {
                  start: function() {
                    var e = [];
                    l && e.push(l.start()), c && e.push(c.start()), z(p, (function(t) {
                      e.push(t.start())
                    }));
                    var t = new o({
                      end: n,
                      cancel: n
                    });
                    return o.all(e, (function(e) {
                      t.complete(e)
                    })), t;

                    function n() {
                      z(e, (function(e) {
                        e.end()
                      }))
                    }
                  }
                }
              }(i.from, i.to, i.classes, i.anchors) : m(i)
            };

            function h(e) {
              return e.replace(/\bng-\S+\b/g, "")
            }

            function g(e, t) {
              return Z(e) && (e = e.split(" ")), Z(t) && (t = t.split(" ")), e.filter((function(e) {
                return -1 === t.indexOf(e)
              })).join(" ")
            }

            function m(t) {
              var n = t.element,
                r = t.options || {};
              t.structural && (r.event = t.event, r.structural = !0, r.applyClassesEarly = !0, "leave" === t.event && (r.onDone = r.domOperation)),
                r.preparationClasses && (r.event = P(r.event, r.preparationClasses));
              var i = e(n, r);
              return i.$$willAnimate ? i : null
            }
          }]
        }],
        oe = ["$animateProvider", function(e) {
          this.$get = ["$injector", "$$AnimateRunner", "$$jqLite", function(t, n, r) {
            var i = S(r);
            return function(e, t, r, a) {
              var s = !1;
              3 === arguments.length && J(r) && (a = r, r = null), a = A(a), r || (r = e.attr("class") || "", a.addClass && (r += " " + a
                .addClass), a.removeClass && (r += " " + a.removeClass));
              var u, l, c, d, p, f = a.addClass,
                h = a.removeClass,
                g = o(r);
              g.length && ("leave" === t ? (d = "leave", c = "afterLeave") : (d = "before" + t.charAt(0).toUpperCase() + t.substr(1), c = t),
                "enter" !== t && "move" !== t && (u = y(e, t, a, g, d)), l = y(e, t, a, g, c));
              if (u || l) return {
                $$willAnimate: !0,
                end: function() {
                  return p ? p.end() : (v(), (p = new n).complete(!0)), p
                },
                start: function() {
                  if (p) return p;
                  var e;
                  p = new n;
                  var t = [];
                  return u && t.push((function(t) {
                    e = u(t)
                  })), t.length ? t.push((function(e) {
                    m(), e(!0)
                  })) : m(), l && t.push((function(t) {
                    e = l(t)
                  })), p.setHost({
                    end: function() {
                      i()
                    },
                    cancel: function() {
                      i(!0)
                    }
                  }), n.chain(t, r), p;

                  function r(e) {
                    v(), p.complete(e)
                  }

                  function i(t) {
                    s || ((e || ne)(t), r(t))
                  }
                }
              };

              function m() {
                a.domOperation(), i(e, a)
              }

              function v() {
                s = !0, m(), T(e, a)
              }

              function $(e, t, r, i, o) {
                var a;
                switch (r) {
                  case "animate":
                    a = [t, i.from, i.to, o];
                    break;
                  case "setClass":
                    a = [t, f, h, o];
                    break;
                  case "addClass":
                    a = [t, f, o];
                    break;
                  case "removeClass":
                    a = [t, h, o];
                    break;
                  default:
                    a = [t, o]
                }
                a.push(i);
                var s = e.apply(e, a);
                if (s)
                  if (X(s.start) && (s = s.start()), s instanceof n) s.done(o);
                  else if (X(s)) return s;
                return ne
              }

              function b(e, t, r, i, o) {
                var a = [];
                return z(i, (function(i) {
                  var s = i[o];
                  s && a.push((function() {
                    var i, o, a = !1,
                      u = function(e) {
                        a || (a = !0, (o || ne)(e), i.complete(!e))
                      };
                    return i = new n({
                      end: function() {
                        u()
                      },
                      cancel: function() {
                        u(!0)
                      }
                    }), o = $(s, e, t, r, (function(e) {
                      u(!1 === e)
                    })), i
                  }))
                })), a
              }

              function y(e, t, r, i, o) {
                var a, s, u = b(e, t, r, i, o);
                0 === u.length && ("beforeSetClass" === o ? (a = b(e, "removeClass", r, i, "beforeRemoveClass"), s = b(e, "addClass", r, i,
                  "beforeAddClass")) : "setClass" === o && (a = b(e, "removeClass", r, i, "removeClass"), s = b(e, "addClass", r, i,
                  "addClass")), a && (u = u.concat(a)), s && (u = u.concat(s)));
                if (0 !== u.length) return function(e) {
                  var t = [];
                  return u.length && z(u, (function(e) {
                      t.push(e())
                    })), t.length ? n.all(t, e) : e(),
                    function(e) {
                      z(t, (function(t) {
                        e ? t.cancel() : t.end()
                      }))
                    }
                }
              }
            };

            function o(n) {
              n = Y(n) ? n : n.split(" ");
              for (var r = [], i = {}, o = 0; o < n.length; o++) {
                var a = n[o],
                  s = e.$$registeredAnimations[a];
                s && !i[a] && (r.push(t.get(s)), i[a] = !0)
              }
              return r
            }
          }]
        }],
        ae = ["$$animationProvider", function(e) {
          e.drivers.push("$$animateJsDriver"), this.$get = ["$$animateJs", "$$AnimateRunner", function(e, t) {
            return function(e) {
              if (e.from && e.to) {
                var r = n(e.from),
                  i = n(e.to);
                if (!r && !i) return;
                return {
                  start: function() {
                    var e = [];
                    r && e.push(r.start()), i && e.push(i.start()), t.all(e, (function(e) {
                      n.complete(e)
                    }));
                    var n = new t({
                      end: o(),
                      cancel: o()
                    });
                    return n;

                    function o() {
                      return function() {
                        z(e, (function(e) {
                          e.end()
                        }))
                      }
                    }
                  }
                }
              }
              return n(e)
            };

            function n(t) {
              var n = t.element,
                r = t.event,
                i = t.options,
                o = t.classes;
              return e(n, r, o, i)
            }
          }]
        }],
        se = "data-ng-animate",
        ue = "$ngAnimatePin",
        le = ["$animateProvider", function(t) {
          var n = this.rules = {
            skip: [],
            cancel: [],
            join: []
          };

          function r(e) {
            return {
              addClass: e.addClass,
              removeClass: e.removeClass,
              from: e.from,
              to: e.to
            }
          }

          function i(e, t) {
            if (e && t) {
              var n = function(e) {
                if (!e) return null;
                var t = e.split(" "),
                  n = Object.create(null);
                return z(t, (function(e) {
                  n[e] = !0
                })), n
              }(t);
              return e.split(" ").some((function(e) {
                return n[e]
              }))
            }
          }

          function o(e, t, r) {
            return n[e].some((function(e) {
              return e(t, r)
            }))
          }

          function u(e, t) {
            var n = (e.addClass || "").length > 0,
              r = (e.removeClass || "").length > 0;
            return t ? n && r : n || r
          }
          n.join.push((function(e, t) {
            return !e.structural && u(e)
          })), n.skip.push((function(e, t) {
            return !e.structural && !u(e)
          })), n.skip.push((function(e, t) {
            return "leave" === t.event && e.structural
          })), n.skip.push((function(e, t) {
            return t.structural && 2 === t.state && !e.structural
          })), n.cancel.push((function(e, t) {
            return t.structural && e.structural
          })), n.cancel.push((function(e, t) {
            return 2 === t.state && e.structural
          })), n.cancel.push((function(e, t) {
            if (t.structural) return !1;
            var n = e.addClass,
              r = e.removeClass,
              o = t.addClass,
              a = t.removeClass;
            return !(ee(n) && ee(r) || ee(o) && ee(a)) && (i(n, a) || i(r, o))
          })), this.$get = ["$$rAF", "$rootScope", "$rootElement", "$document", "$$Map", "$$animation", "$$AnimateRunner", "$templateRequest",
            "$$jqLite", "$$forceReflow", "$$isDocumentHidden",
            function(n, i, c, d, p, f, h, g, m, v, $) {
              var y = new p,
                k = new p,
                D = null;

              function O(e) {
                k.delete(e.target)
              }
              var N = i.$watch((function() {
                  return 0 === g.totalPendingRequests
                }), (function(e) {
                  e && (N(), i.$$postDigest((function() {
                    i.$$postDigest((function() {
                      null === D && (D = !0)
                    }))
                  })))
                })),
                q = Object.create(null),
                I = t.customFilter(),
                L = t.classNameFilter(),
                R = function() {
                  return !0
                },
                j = I || R,
                _ = L ? function(e, t) {
                  var n = [e.getAttribute("class"), t.addClass, t.removeClass].join(" ");
                  return L.test(n)
                } : R,
                V = S(m);

              function U(e, t) {
                return M(e, t, {})
              }
              var H = e.Node.prototype.contains || function(e) {
                return this === e || !!(16 & this.compareDocumentPosition(e))
              };

              function F(e, t, n) {
                var r = C(t);
                return e.filter((function(e) {
                  return !(e.node === r && (!n || e.callback === n))
                }))
              }

              function B(e, t) {
                "close" !== e || t.parentNode || X.off(t)
              }
              var X = {
                on: function(e, t, n) {
                  var r = C(t);
                  q[e] = q[e] || [], q[e].push({
                    node: r,
                    callback: n
                  }), te(t).on("$destroy", (function() {
                    y.get(r) || X.off(e, t, n)
                  }))
                },
                off: function(e, t, n) {
                  if (1 !== arguments.length || Z(arguments[0])) {
                    var r = q[e];
                    r && (q[e] = 1 === arguments.length ? null : F(r, t, n))
                  } else
                    for (var i in t = arguments[0], q) q[i] = F(q[i], t)
                },
                pin: function(e, t) {
                  b(Q(e), "element", "not an element"), b(Q(t), "parentElement", "not an element"), e.data(ue, t)
                },
                push: function(e, t, p, g) {
                  return (p = p || {}).domOperation = g,
                    function(e, t, p) {
                      var g = G(p),
                        m = x(e),
                        v = E(m),
                        b = v && v.parentNode;
                      g = A(g);
                      var C = new h,
                        S = (O = !1, function(e) {
                          O ? e() : i.$$postDigest((function() {
                            O = !0, e()
                          }))
                        });
                      var O;
                      Y(g.addClass) && (g.addClass = g.addClass.join(" "));
                      g.addClass && !Z(g.addClass) && (g.addClass = null);
                      Y(g.removeClass) && (g.removeClass = g.removeClass.join(" "));
                      g.removeClass && !Z(g.removeClass) && (g.removeClass = null);
                      g.from && !J(g.from) && (g.from = null);
                      g.to && !J(g.to) && (g.to = null);
                      if (!(D && v && j(v, t, p) && _(v, g))) return oe(), C;
                      var N = ["enter", "move", "leave"].indexOf(t) >= 0,
                        I = $(),
                        L = I || k.get(v),
                        R = !L && y.get(v) || {},
                        F = !!R.state;
                      L || F && 1 === R.state || (L = ! function(e, t, n) {
                        var r, i = d[0].body,
                          o = E(c),
                          a = e === i || "HTML" === e.nodeName,
                          s = e === o,
                          u = !1,
                          p = k.get(e),
                          f = te.data(e, ue);
                        f && (t = E(f));
                        for (; t && (s || (s = t === o), 1 === t.nodeType);) {
                          var h = y.get(t) || {};
                          if (!u) {
                            var g = k.get(t);
                            if (!0 === g && !1 !== p) {
                              p = !0;
                              break
                            }!1 === g && (p = !1), u = h.structural
                          }
                          if (ee(r) || !0 === r) {
                            var m = te.data(t, l);
                            K(m) && (r = m)
                          }
                          if (u && !1 === r) break;
                          if (a || (a = t === i), a && s) break;
                          t = s || !(f = te.data(t, ue)) ? t.parentNode : E(f)
                        }
                        return (!u || r) && !0 !== p && s && a
                      }(v, b));
                      if (L) return I && ie(C, t, "start", r(g)), oe(), I && ie(C, t, "close", r(g)), C;
                      N && function(e) {
                        var t = e.querySelectorAll("[data-ng-animate]");
                        z(t, (function(e) {
                          var t = parseInt(e.getAttribute(se), 10),
                            n = y.get(e);
                          if (n) switch (t) {
                            case 2:
                              n.runner.end();
                            case 1:
                              y.delete(e)
                          }
                        }))
                      }(v);
                      var W = {
                        structural: N,
                        element: m,
                        event: t,
                        addClass: g.addClass,
                        removeClass: g.removeClass,
                        close: oe,
                        options: g,
                        runner: C
                      };
                      if (F) {
                        if (o("skip", W, R)) return 2 === R.state ? (oe(), C) : (M(m, R, W), R.runner);
                        if (o("cancel", W, R))
                          if (2 === R.state) R.runner.end();
                          else {
                            if (!R.structural) return M(m, R, W), R.runner;
                            R.close()
                          }
                        else if (o("join", W, R)) {
                          if (2 !== R.state) return function(e, t, n, r) {
                            var i = "";
                            n && (i = w(n, "ng-", !0)), r.addClass && (i = P(i, w(r.addClass, a))), r.removeClass && (i = P(i, w(r
                              .removeClass, s))), i.length && (r.preparationClasses = i, t.addClass(i))
                          }(0, m, N ? t : null, g), t = W.event = R.event, g = M(m, R, W), R.runner;
                          U(m, W)
                        }
                      } else U(m, W);
                      var Q = W.structural;
                      Q || (Q = "animate" === W.event && Object.keys(W.options.to || {}).length > 0 || u(W));
                      if (!Q) return oe(), ne(v), C;
                      var X = (R.counter || 0) + 1;
                      return W.counter = X, re(v, 1, W), i.$$postDigest((function() {
                        m = x(e);
                        var n = y.get(v),
                          i = !n;
                        n = n || {};
                        var o = (m.parent() || []).length > 0 && ("animate" === n.event || n.structural || u(n));
                        if (i || n.counter !== X || !o) return i && (V(m, g), T(m, g)), (i || N && n.event !== t) && (g.domOperation(), C
                          .end()), void(o || ne(v));
                        t = !n.structural && u(n, !0) ? "setClass" : n.event, re(v, 2);
                        var a = f(m, t, n.options);
                        C.setHost(a), ie(C, t, "start", r(g)), a.done((function(e) {
                          oe(!e);
                          var n = y.get(v);
                          n && n.counter === X && ne(v), ie(C, t, "close", r(g))
                        }))
                      })), C;

                      function ie(e, t, r, i) {
                        S((function() {
                          var e = function(e, t, n) {
                            var r = [],
                              i = q[n];
                            return i && z(i, (function(i) {
                              (H.call(i.node, t) || "leave" === n && H.call(i.node, e)) && r.push(i.callback)
                            })), r
                          }(b, v, t);
                          e.length ? n((function() {
                            z(e, (function(e) {
                              e(m, r, i)
                            })), B(r, v)
                          })) : B(r, v)
                        })), e.progress(t, r, i)
                      }

                      function oe(e) {
                        ! function(e, t) {
                          t.preparationClasses && (e.removeClass(t.preparationClasses), t.preparationClasses = null), t.activeClasses && (e
                            .removeClass(t.activeClasses), t.activeClasses = null)
                        }(m, g), V(m, g), T(m, g), g.domOperation(), C.complete(!e)
                      }
                    }(e, t, p)
                },
                enabled: function(e, t) {
                  var n = arguments.length;
                  if (0 === n) t = !!D;
                  else {
                    var r = Q(e);
                    if (r) {
                      var i = E(e);
                      1 === n ? t = !k.get(i) : (k.has(i) || te(e).on("$destroy", O), k.set(i, !t))
                    } else t = D = !!e
                  }
                  return t
                }
              };
              return X;

              function ne(e) {
                e.removeAttribute(se), y.delete(e)
              }

              function re(e, t, n) {
                (n = n || {}).state = t, e.setAttribute(se, t);
                var r = y.get(e),
                  i = r ? W(r, n) : n;
                y.set(e, i)
              }
            }
          ]
        }],
        ce = ["$animateProvider", function(e) {
          var t = "ng-animate-ref",
            n = this.drivers = [],
            r = "$$animationRunner",
            i = "$$animatePrepareClasses";

          function o(e) {
            return e.data(r)
          }
          this.$get = ["$$jqLite", "$rootScope", "$injector", "$$AnimateRunner", "$$Map", "$$rAFScheduler", "$$animateCache", function(e, a, s, l, c, d,
            p) {
            var f = [],
              h = S(e);
            return function(g, m, v) {
              v = A(v);
              var $ = ["enter", "move", "leave"].indexOf(m) >= 0,
                b = new l({
                  end: function() {
                    k()
                  },
                  cancel: function() {
                    k(!0)
                  }
                });
              if (!n.length) return k(), b;
              var w = y(g.attr("class"), y(v.addClass, v.removeClass)),
                x = v.tempClasses;
              return x && (w += " " + x, v.tempClasses = null), $ && g.data(i, "ng-" + m + "-prepare"),
                function(e, t) {
                  e.data(r, t)
                }(g, b), f.push({
                  element: g,
                  classes: w,
                  event: m,
                  structural: $,
                  options: v,
                  beforeStart: function() {
                    x = (x ? x + " " : "") + u, e.addClass(g, x);
                    var t = g.data(i);
                    t && (e.removeClass(g, t), t = null)
                  },
                  close: k
                }), g.on("$destroy", S), f.length > 1 || a.$$postDigest((function() {
                  var r = [];
                  z(f, (function(e) {
                    o(e.element) ? r.push(e) : e.close()
                  })), f.length = 0;
                  var a = function(e) {
                      var n = [],
                        r = {};
                      z(e, (function(e, i) {
                        var o = E(e.element),
                          a = e.event,
                          s = ["enter", "move"].indexOf(a) >= 0,
                          u = e.structural ? function(e) {
                            var n = "[ng-animate-ref]",
                              r = e.hasAttribute(t) ? [e] : e.querySelectorAll(n),
                              i = [];
                            return z(r, (function(e) {
                              var n = e.getAttribute(t);
                              n && n.length && i.push(e)
                            })), i
                          }(o) : [];
                        if (u.length) {
                          var l = s ? "to" : "from";
                          z(u, (function(e) {
                            var n = e.getAttribute(t);
                            r[n] = r[n] || {}, r[n][l] = {
                              animationID: i,
                              element: te(e)
                            }
                          }))
                        } else n.push(e)
                      }));
                      var i = {},
                        o = {};
                      return z(r, (function(t, r) {
                        var a = t.from,
                          s = t.to;
                        if (a && s) {
                          var u = e[a.animationID],
                            l = e[s.animationID],
                            c = a.animationID.toString();
                          if (!o[c]) {
                            var d = o[c] = {
                              structural: !0,
                              beforeStart: function() {
                                u.beforeStart(), l.beforeStart()
                              },
                              close: function() {
                                u.close(), l.close()
                              },
                              classes: C(u.classes, l.classes),
                              from: u,
                              to: l,
                              anchors: []
                            };
                            d.classes.length ? n.push(d) : (n.push(u), n.push(l))
                          }
                          o[c].anchors.push({
                            out: a.element,
                            in: s.element
                          })
                        } else {
                          var p = a ? a.animationID : s.animationID,
                            f = p.toString();
                          i[f] || (i[f] = !0, n.push(e[p]))
                        }
                      })), n
                    }(r),
                    l = [];
                  z(a, (function(e) {
                    var t = e.from ? e.from.element : e.element,
                      r = v.addClass;
                    r = (r ? r + " " : "") + u;
                    var i = p.cacheKey(t[0], e.event, r, v.removeClass);
                    l.push({
                      element: t,
                      domNode: E(t),
                      fn: function() {
                        var t, r = e.close;
                        if (p.containsCachedAnimationWithoutDuration(i)) r();
                        else {
                          if (e.beforeStart(), o(e.anchors ? e.from.element || e.to.element : e.element)) {
                            var a = function(e) {
                              for (var t = n.length - 1; t >= 0; t--) {
                                var r = n[t],
                                  i = s.get(r)(e);
                                if (i) return i
                              }
                            }(e);
                            a && (t = a.start)
                          }
                          if (t) {
                            var u = t();
                            u.done((function(e) {
                                r(!e)
                              })),
                              function(e, t) {
                                e.from && e.to ? (n(e.from.element), n(e.to.element)) : n(e.element);

                                function n(e) {
                                  var n = o(e);
                                  n && n.setHost(t)
                                }
                              }(e, u)
                          } else r()
                        }
                      }
                    })
                  }));
                  for (var h = function(e) {
                      var t, n = {
                          children: []
                        },
                        r = new c;
                      for (t = 0; t < e.length; t++) {
                        var i = e[t];
                        r.set(i.domNode, e[t] = {
                          domNode: i.domNode,
                          element: i.element,
                          fn: i.fn,
                          children: []
                        })
                      }
                      for (t = 0; t < e.length; t++) o(e[t]);
                      return function(e) {
                        var t, n = [],
                          r = [];
                        for (t = 0; t < e.children.length; t++) r.push(e.children[t]);
                        var i = r.length,
                          o = 0,
                          a = [];
                        for (t = 0; t < r.length; t++) {
                          var s = r[t];
                          i <= 0 && (i = o, o = 0, n.push(a), a = []), a.push(s), s.children.forEach((function(e) {
                            o++, r.push(e)
                          })), i--
                        }
                        return a.length && n.push(a), n
                      }(n);

                      function o(e) {
                        if (e.processed) return e;
                        e.processed = !0;
                        var t, i = e.domNode,
                          a = i.parentNode;
                        for (r.set(i, e); a;) {
                          if (t = r.get(a)) {
                            t.processed || (t = o(t));
                            break
                          }
                          a = a.parentNode
                        }
                        return (t || n).children.push(e), e
                      }
                    }(l), g = 0; g < h.length; g++)
                    for (var m = h[g], $ = 0; $ < m.length; $++) {
                      var b = m[$],
                        y = b.element;
                      if (h[g][$] = b.fn, 0 !== g) {
                        var w = y.data(i);
                        w && e.addClass(y, w)
                      } else y.removeData(i)
                    }
                  d(h)
                })), b;

              function C(e, t) {
                e = e.split(" "), t = t.split(" ");
                for (var n = [], r = 0; r < e.length; r++) {
                  var i = e[r];
                  if ("ng-" !== i.substring(0, 3))
                    for (var o = 0; o < t.length; o++)
                      if (i === t[o]) {
                        n.push(i);
                        break
                      }
                }
                return n.join(" ")
              }

              function S() {
                var e = o(g);
                !e || "leave" === m && v.$$domOperationFired || e.end()
              }

              function k(t) {
                g.off("$destroy", S),
                  function(e) {
                    e.removeData(r)
                  }(g), h(g, v), T(g, v), v.domOperation(), x && e.removeClass(g, x), b.complete(!t)
              }
            }
          }]
        }];
      t.module("ngAnimate", [], (function() {
        ne = t.noop, G = t.copy, W = t.extend, te = t.element, z = t.forEach, Y = t.isArray, Z = t.isString, J = t.isObject, ee = t.isUndefined, K =
          t.isDefined, X = t.isFunction, Q = t.isElement
      })).info({
        angularVersion: "1.8.2"
      }).directive("ngAnimateSwap", ["$animate", function(e) {
        return {
          restrict: "A",
          transclude: "element",
          terminal: !0,
          priority: 550,
          link: function(t, n, r, i, o) {
            var a, s;
            t.$watchCollection(r.ngAnimateSwap || r.for, (function(t) {
              a && e.leave(a), s && (s.$destroy(), s = null), (t || 0 === t) && o((function(t, r) {
                a = t, s = r, e.enter(t, null, n)
              }))
            }))
          }
        }
      }]).directive("ngAnimateChildren", I).factory("$$rAFScheduler", ["$$rAF", function(e) {
        var t, n;

        function r(e) {
          t = t.concat(e), i()
        }
        return t = r.queue = [], r.waitUntilQuiet = function(t) {
          n && n(), n = e((function() {
            n = null, t(), i()
          }))
        }, r;

        function i() {
          if (t.length) {
            for (var r = t.shift(), o = 0; o < r.length; o++) r[o]();
            n || e((function() {
              n || i()
            }))
          }
        }
      }]).provider("$$animateQueue", le).provider("$$animateCache", (function() {
        var e = "$$ngAnimateParentKey",
          t = 0,
          n = Object.create(null);
        this.$get = [function() {
          return {
            cacheKey: function(n, r, i, o) {
              var a = n.parentNode,
                s = [a[e] || (a[e] = ++t), r, n.getAttribute("class")];
              return i && s.push(i), o && s.push(o), s.join(" ")
            },
            containsCachedAnimationWithoutDuration: function(e) {
              var t = n[e];
              return t && !t.isValid || !1
            },
            flush: function() {
              n = Object.create(null)
            },
            count: function(e) {
              var t = n[e];
              return t ? t.total : 0
            },
            get: function(e) {
              var t = n[e];
              return t && t.value
            },
            put: function(e, t, r) {
              n[e] ? (n[e].total++, n[e].value = t) : n[e] = {
                total: 1,
                value: t,
                isValid: r
              }
            }
          }
        }]
      })).provider("$$animation", ce).provider("$animateCss", re).provider("$$animateCssDriver", ie).provider("$$animateJs", oe).provider(
        "$$animateJsDriver", ae)
    }(window, window.angular)
  }, {}],
  24: [function(e, t, n) {
    e("./angular-animate"), t.exports = "ngAnimate"
  }, {
    "./angular-animate": 23
  }],
  25: [function(e, t, n) {
    /**
     * @license AngularJS v1.8.2
     * (c) 2010-2020 Google LLC. http://angularjs.org
     * License: MIT
     */
    ! function(e, t) {
      "use strict";

      function n(e, t) {
        var n = [],
          r = e.replace(/([().])/g, "\\$1").replace(/(\/)?:(\w+)(\*\?|[?*])?/g, (function(e, t, r, i) {
            var o = "?" === i || "*?" === i,
              a = "*" === i || "*?" === i;
            return n.push({
              name: r,
              optional: o
            }), t = t || "", (o ? "(?:" + t : t + "(?:") + (a ? "(.+?)" : "([^/]+)") + (o ? "?)?" : ")")
          })).replace(/([/$*])/g, "\\$1");
        return t.ignoreTrailingSlashes && (r = r.replace(/\/+$/, "") + "/*"), {
          keys: n,
          regexp: new RegExp("^" + r + "(?:[?#]|$)", t.caseInsensitiveMatch ? "i" : "")
        }
      }
      var r, i, o, a, s, u = t.module("ngRoute", []).info({
          angularVersion: "1.8.2"
        }).provider("$route", (function() {
          function e(e, n) {
            return t.extend(Object.create(e), n)
          }
          r = t.isArray, i = t.isObject, o = t.isDefined, a = t.noop;
          var u = {};
          this.when = function(e, o) {
            var a = function(e, t) {
              if (r(e)) {
                t = t || [];
                for (var n = 0, o = e.length; n < o; n++) t[n] = e[n]
              } else if (i(e))
                for (var a in t = t || {}, e) "$" === a.charAt(0) && "$" === a.charAt(1) || (t[a] = e[a]);
              return t || e
            }(o);
            if (t.isUndefined(a.reloadOnUrl) && (a.reloadOnUrl = !0), t.isUndefined(a.reloadOnSearch) && (a.reloadOnSearch = !0), t.isUndefined(a
                .caseInsensitiveMatch) && (a.caseInsensitiveMatch = this.caseInsensitiveMatch), u[e] = t.extend(a, {
                originalPath: e
              }, e && n(e, a)), e) {
              var s = "/" === e[e.length - 1] ? e.substr(0, e.length - 1) : e + "/";
              u[s] = t.extend({
                originalPath: e,
                redirectTo: e
              }, n(s, a))
            }
            return this
          }, this.caseInsensitiveMatch = !1, this.otherwise = function(e) {
            return "string" == typeof e && (e = {
              redirectTo: e
            }), this.when(null, e), this
          }, s = !0, this.eagerInstantiationEnabled = function(e) {
            return o(e) ? (s = e, this) : s
          }, this.$get = ["$rootScope", "$location", "$routeParams", "$q", "$injector", "$templateRequest", "$sce", "$browser", function(n, r, i, o,
            s, c, d, p) {
            var f, h, g = !1,
              m = {
                routes: u,
                reload: function() {
                  g = !0;
                  var e = {
                    defaultPrevented: !1,
                    preventDefault: function() {
                      this.defaultPrevented = !0, g = !1
                    }
                  };
                  n.$evalAsync((function() {
                    v(e), e.defaultPrevented || $()
                  }))
                },
                updateParams: function(e) {
                  if (!this.current || !this.current.$$route) throw l("norout", "Tried updating route with no current route");
                  e = t.extend({}, this.current.params, e), r.path(x(this.current.$$route.originalPath, e)), r.search(e)
                }
              };
            return n.$on("$locationChangeStart", v), n.$on("$locationChangeSuccess", $), m;

            function v(i) {
              var o, a, s, l, c = m.current;
              t.forEach(u, (function(n, i) {
                !a && (o = function(e, t) {
                  var n = t.keys,
                    r = {};
                  if (!t.regexp) return null;
                  var i = t.regexp.exec(e);
                  if (!i) return null;
                  for (var o = 1, a = i.length; o < a; ++o) {
                    var s = n[o - 1],
                      u = i[o];
                    s && u && (r[s.name] = u)
                  }
                  return r
                }(r.path(), n)) && ((a = e(n, {
                  params: t.extend({}, r.search(), o),
                  pathParams: o
                })).$$route = n)
              })), f = a || u.null && e(u.null, {
                params: {},
                pathParams: {}
              }), s = f, l = c, (h = !g && s && l && s.$$route === l.$$route && (!s.reloadOnUrl || !s.reloadOnSearch && t.equals(s.pathParams, l
                .pathParams))) || !c && !f || n.$broadcast("$routeChangeStart", f, c).defaultPrevented && i && i.preventDefault()
            }

            function $() {
              var e = m.current,
                r = f;
              if (h) e.params = r.params, t.copy(e.params, i), n.$broadcast("$routeUpdate", e);
              else if (r || e) {
                g = !1, m.current = r;
                var s = o.resolve(r);
                p.$$incOutstandingRequestCount("$route"), s.then(b).then(y).then((function(o) {
                  return o && s.then(w).then((function(o) {
                    r === m.current && (r && (r.locals = o, t.copy(r.params, i)), n.$broadcast("$routeChangeSuccess", r, e))
                  }))
                })).catch((function(t) {
                  r === m.current && n.$broadcast("$routeChangeError", r, e, t)
                })).finally((function() {
                  p.$$completeOutstandingRequest(a, "$route")
                }))
              }
            }

            function b(e) {
              var n = {
                route: e,
                hasRedirection: !1
              };
              if (e)
                if (e.redirectTo)
                  if (t.isString(e.redirectTo)) n.path = x(e.redirectTo, e.params), n.search = e.params, n.hasRedirection = !0;
                  else {
                    var i = r.path(),
                      a = r.search(),
                      u = e.redirectTo(e.pathParams, i, a);
                    t.isDefined(u) && (n.url = u, n.hasRedirection = !0)
                  }
              else if (e.resolveRedirectTo) return o.resolve(s.invoke(e.resolveRedirectTo)).then((function(e) {
                return t.isDefined(e) && (n.url = e, n.hasRedirection = !0), n
              }));
              return n
            }

            function y(e) {
              var t = !0;
              if (e.route !== m.current) t = !1;
              else if (e.hasRedirection) {
                var n = r.url(),
                  i = e.url;
                i ? r.url(i).replace() : i = r.path(e.path).search(e.search).replace().url(), i !== n && (t = !1)
              }
              return t
            }

            function w(e) {
              if (e) {
                var n = t.extend({}, e.resolve);
                t.forEach(n, (function(e, r) {
                  n[r] = t.isString(e) ? s.get(e) : s.invoke(e, null, null, r)
                }));
                var r = function(e) {
                  var n, r;
                  t.isDefined(n = e.template) ? t.isFunction(n) && (n = n(e.params)) : t.isDefined(r = e.templateUrl) && (t.isFunction(r) && (
                    r = r(e.params)), t.isDefined(r) && (e.loadedTemplateUrl = d.valueOf(r), n = c(r)));
                  return n
                }(e);
                return t.isDefined(r) && (n.$template = r), o.all(n)
              }
            }

            function x(e, n) {
              var r = [];
              return t.forEach((e || "").split(":"), (function(e, t) {
                if (0 === t) r.push(e);
                else {
                  var i = e.match(/(\w+)(?:[?*])?(.*)/),
                    o = i[1];
                  r.push(n[o]), r.push(i[2] || ""), delete n[o]
                }
              })), r.join("")
            }
          }]
        })).run(c),
        l = t.$$minErr("ngRoute");

      function c(e) {
        s && e.get("$route")
      }

      function d(e, n, r) {
        return {
          restrict: "ECA",
          terminal: !0,
          priority: 400,
          transclude: "element",
          link: function(i, o, a, s, u) {
            var l, c, d, p = a.autoscroll,
              f = a.onload || "";

            function h() {
              d && (r.cancel(d), d = null), l && (l.$destroy(), l = null), c && ((d = r.leave(c)).done((function(e) {
                !1 !== e && (d = null)
              })), c = null)
            }

            function g() {
              var a = e.current && e.current.locals,
                s = a && a.$template;
              if (t.isDefined(s)) {
                var d = i.$new(),
                  g = e.current,
                  m = u(d, (function(e) {
                    r.enter(e, null, c || o).done((function(e) {
                      !1 === e || !t.isDefined(p) || p && !i.$eval(p) || n()
                    })), h()
                  }));
                c = m, (l = g.scope = d).$emit("$viewContentLoaded"), l.$eval(f)
              } else h()
            }
            i.$on("$routeChangeSuccess", g), g()
          }
        }
      }

      function p(e, t, n) {
        return {
          restrict: "ECA",
          priority: -400,
          link: function(r, i) {
            var o = n.current,
              a = o.locals;
            i.html(a.$template);
            var s = e(i.contents());
            if (o.controller) {
              a.$scope = r;
              var u = t(o.controller, a);
              o.controllerAs && (r[o.controllerAs] = u), i.data("$ngControllerController", u), i.children().data("$ngControllerController", u)
            }
            r[o.resolveAs || "$resolve"] = a, s(r)
          }
        }
      }
      c.$inject = ["$injector"], u.provider("$routeParams", (function() {
        this.$get = function() {
          return {}
        }
      })), u.directive("ngView", d), u.directive("ngView", p), d.$inject = ["$route", "$anchorScroll", "$animate"], p.$inject = ["$compile",
        "$controller", "$route"
      ]
    }(window, window.angular)
  }, {}],
  26: [function(e, t, n) {
    e("./angular-route"), t.exports = "ngRoute"
  }, {
    "./angular-route": 25
  }],
  27: [function(e, t, n) {
    /**
     * @license AngularJS v1.8.2
     * (c) 2010-2020 Google LLC. http://angularjs.org
     * License: MIT
     */
    ! function(e, t) {
      "use strict";
      var n, r, i, o, a, s, u, l, c, d, p = t.$$minErr("$sanitize");
      t.module("ngSanitize", []).provider("$sanitize", (function() {
        var f = !1,
          h = !1;
        this.$get = ["$$sanitizeUri", function(e) {
          return f = !0, h && r(A, C),
            function(t) {
              var n = [];
              return c(t, d(n, (function(t, n) {
                return !/^unsafe:/.test(e(t, n))
              }))), n.join("")
            }
        }], this.enableSvg = function(e) {
          return a(e) ? (h = e, this) : h
        }, this.addValidElements = function(e) {
          return f || (o(e) && (e = {
            htmlElements: e
          }), N(C, e.svgElements), N(v, e.htmlVoidElements), N(A, e.htmlVoidElements), N(A, e.htmlElements)), this
        }, this.addValidAttrs = function(e) {
          return f || r(M, O(e, !0)), this
        }, n = t.bind, r = t.extend, i = t.forEach, o = t.isArray, a = t.isDefined, s = t.$$lowercase, u = t.noop, c = function(e, t) {
          null == e ? e = "" : "string" != typeof e && (e = "" + e);
          var n = P(e);
          if (!n) return "";
          var r = 5;
          do {
            if (0 === r) throw p("uinput", "Failed to sanitize html because the input is unstable");
            r--, e = n.innerHTML, n = P(e)
          } while (e !== n.innerHTML);
          var i = n.firstChild;
          for (; i;) {
            switch (i.nodeType) {
              case 1:
                t.start(i.nodeName.toLowerCase(), q(i.attributes));
                break;
              case 3:
                t.chars(i.textContent)
            }
            var o;
            if (!((o = i.firstChild) || (1 === i.nodeType && t.end(i.nodeName.toLowerCase()), o = R("nextSibling", i))))
              for (; null == o && (i = R("parentNode", i)) !== n;) o = R("nextSibling", i), 1 === i.nodeType && t.end(i.nodeName.toLowerCase());
            i = o
          }
          for (; i = n.firstChild;) n.removeChild(i)
        }, d = function(e, t) {
          var r = !1,
            o = n(e, e.push);
          return {
            start: function(e, n) {
              e = s(e), !r && S[e] && (r = e), r || !0 !== A[e] || (o("<"), o(e), i(n, (function(n, r) {
                var i = s(r),
                  a = "img" === e && "src" === i || "background" === i;
                !0 !== M[i] || !0 === T[i] && !t(n, a) || (o(" "), o(r), o('="'), o(I(n)), o('"'))
              })), o(">"))
            },
            end: function(e) {
              e = s(e), r || !0 !== A[e] || !0 === v[e] || (o("</"), o(e), o(">")), e == r && (r = !1)
            },
            chars: function(e) {
              r || o(I(e))
            }
          }
        }, l = e.Node.prototype.contains || function(e) {
          return !!(16 & this.compareDocumentPosition(e))
        };
        var g = /[\uD800-\uDBFF][\uDC00-\uDFFF]/g,
          m = /([^#-~ |!])/g,
          v = E("area,br,col,hr,img,wbr"),
          $ = E("colgroup,dd,dt,li,p,tbody,td,tfoot,th,thead,tr"),
          b = E("rp,rt"),
          y = r({}, b, $),
          w = r({}, $, E(
            "address,article,aside,blockquote,caption,center,del,dir,div,dl,figure,figcaption,footer,h1,h2,h3,h4,h5,h6,header,hgroup,hr,ins,map,menu,nav,ol,pre,section,table,ul"
            )),
          x = r({}, b, E(
            "a,abbr,acronym,b,bdi,bdo,big,br,cite,code,del,dfn,em,font,i,img,ins,kbd,label,map,mark,q,ruby,rp,rt,s,samp,small,span,strike,strong,sub,sup,time,tt,u,var"
            )),
          C = E(
            "circle,defs,desc,ellipse,font-face,font-face-name,font-face-src,g,glyph,hkern,image,linearGradient,line,marker,metadata,missing-glyph,mpath,path,polygon,polyline,radialGradient,rect,stop,svg,switch,text,title,tspan"
            ),
          S = E("script,style"),
          A = r({}, v, w, x, y),
          T = E("background,cite,href,longdesc,src,xlink:href,xml:base"),
          k = E(
            "abbr,align,alt,axis,bgcolor,border,cellpadding,cellspacing,class,clear,color,cols,colspan,compact,coords,dir,face,headers,height,hreflang,hspace,ismap,lang,language,nohref,nowrap,rel,rev,rows,rowspan,rules,scope,scrolling,shape,size,span,start,summary,tabindex,target,title,type,valign,value,vspace,width"
            ),
          D = E(
            "accent-height,accumulate,additive,alphabetic,arabic-form,ascent,baseProfile,bbox,begin,by,calcMode,cap-height,class,color,color-rendering,content,cx,cy,d,dx,dy,descent,display,dur,end,fill,fill-rule,font-family,font-size,font-stretch,font-style,font-variant,font-weight,from,fx,fy,g1,g2,glyph-name,gradientUnits,hanging,height,horiz-adv-x,horiz-origin-x,ideographic,k,keyPoints,keySplines,keyTimes,lang,marker-end,marker-mid,marker-start,markerHeight,markerUnits,markerWidth,mathematical,max,min,offset,opacity,orient,origin,overline-position,overline-thickness,panose-1,path,pathLength,points,preserveAspectRatio,r,refX,refY,repeatCount,repeatDur,requiredExtensions,requiredFeatures,restart,rotate,rx,ry,slope,stemh,stemv,stop-color,stop-opacity,strikethrough-position,strikethrough-thickness,stroke,stroke-dasharray,stroke-dashoffset,stroke-linecap,stroke-linejoin,stroke-miterlimit,stroke-opacity,stroke-width,systemLanguage,target,text-anchor,to,transform,type,u1,u2,underline-position,underline-thickness,unicode,unicode-range,units-per-em,values,version,viewBox,visibility,width,widths,x,x-height,x1,x2,xlink:actuate,xlink:arcrole,xlink:role,xlink:show,xlink:title,xlink:type,xml:base,xml:lang,xml:space,xmlns,xmlns:xlink,y,y1,y2,zoomAndPan",
            !0),
          M = r({}, T, D, k);

        function E(e, t) {
          return O(e.split(","), t)
        }

        function O(e, t) {
          var n, r = {};
          for (n = 0; n < e.length; n++) r[t ? s(e[n]) : e[n]] = !0;
          return r
        }

        function N(e, t) {
          t && t.length && r(e, O(t))
        }
        var P = function(e, t) {
          if (function() {
              try {
                return !!i("")
              } catch (e) {
                return !1
              }
            }()) return i;
          if (!t || !t.implementation) throw p("noinert", "Can't create an inert html document");
          var n = t.implementation.createHTMLDocument("inert"),
            r = (n.documentElement || n.getDocumentElement()).querySelector("body");
          return function(e) {
            r.innerHTML = e, t.documentMode && L(r);
            return r
          };

          function i(t) {
            t = "<remove></remove>" + t;
            try {
              var n = (new e.DOMParser).parseFromString(t, "text/html").body;
              return n.firstChild.remove(), n
            } catch (e) {
              return
            }
          }
        }(e, e.document);

        function q(e) {
          for (var t = {}, n = 0, r = e.length; n < r; n++) {
            var i = e[n];
            t[i.name] = i.value
          }
          return t
        }

        function I(e) {
          return e.replace(/&/g, "&amp;").replace(g, (function(e) {
            return "&#" + (1024 * (e.charCodeAt(0) - 55296) + (e.charCodeAt(1) - 56320) + 65536) + ";"
          })).replace(m, (function(e) {
            return "&#" + e.charCodeAt(0) + ";"
          })).replace(/</g, "&lt;").replace(/>/g, "&gt;")
        }

        function L(t) {
          for (; t;) {
            if (t.nodeType === e.Node.ELEMENT_NODE)
              for (var n = t.attributes, r = 0, i = n.length; r < i; r++) {
                var o = n[r],
                  a = o.name.toLowerCase();
                "xmlns:ns1" !== a && 0 !== a.lastIndexOf("ns1:", 0) || (t.removeAttributeNode(o), r--, i--)
              }
            var s = t.firstChild;
            s && L(s), t = R("nextSibling", t)
          }
        }

        function R(e, t) {
          var n = t[e];
          if (n && l.call(t, n)) throw p("elclob", "Failed to sanitize html because the element is clobbered: {0}", t.outerHTML || t.outerText);
          return n
        }
      })).info({
        angularVersion: "1.8.2"
      }), t.module("ngSanitize").filter("linky", ["$sanitize", function(e) {
        var n = /((s?ftp|https?):\/\/|(www\.)|(mailto:)?[A-Za-z0-9._%+-]+@)\S*[^\s.;,(){}<>"\u201d\u2019]/i,
          r = /^mailto:/i,
          i = t.$$minErr("linky"),
          o = t.isDefined,
          a = t.isFunction,
          s = t.isObject,
          l = t.isString;
        return function(t, c, p) {
          if (null == t || "" === t) return t;
          if (!l(t)) throw i("notstring", "Expected string but received: {0}", t);
          for (var f, h, g, m = a(p) ? p : s(p) ? function() {
              return p
            } : function() {
              return {}
            }, v = t, $ = []; f = v.match(n);) h = f[0], f[2] || f[4] || (h = (f[3] ? "http://" : "mailto:") + h), g = f.index, b(v.substr(0, g)),
            y(h, f[0].replace(r, "")), v = v.substring(g + f[0].length);
          return b(v), e($.join(""));

          function b(e) {
            var t, n;
            e && $.push((t = e, d(n = [], u).chars(t), n.join("")))
          }

          function y(e, t) {
            var n, r = m(e);
            for (n in $.push("<a "), r) $.push(n + '="' + r[n] + '" ');
            o(c) && !("target" in r) && $.push('target="', c, '" '), $.push('href="', e.replace(/"/g, "&quot;"), '">'), b(t), $.push("</a>")
          }
        }
      }])
    }(window, window.angular)
  }, {}],
  28: [function(e, t, n) {
    e("./angular-sanitize"), t.exports = "ngSanitize"
  }, {
    "./angular-sanitize": 27
  }],
  29: [function(e, t, n) {
    angular.module("ui.bootstrap", ["ui.bootstrap.tpls", "ui.bootstrap.collapse", "ui.bootstrap.tabindex", "ui.bootstrap.accordion", "ui.bootstrap.alert",
      "ui.bootstrap.buttons", "ui.bootstrap.carousel", "ui.bootstrap.dateparser", "ui.bootstrap.isClass", "ui.bootstrap.datepicker",
      "ui.bootstrap.position", "ui.bootstrap.datepickerPopup", "ui.bootstrap.debounce", "ui.bootstrap.multiMap", "ui.bootstrap.dropdown",
      "ui.bootstrap.stackedMap", "ui.bootstrap.modal", "ui.bootstrap.paging", "ui.bootstrap.pager", "ui.bootstrap.pagination", "ui.bootstrap.tooltip",
      "ui.bootstrap.popover", "ui.bootstrap.progressbar", "ui.bootstrap.rating", "ui.bootstrap.tabs", "ui.bootstrap.timepicker",
      "ui.bootstrap.typeahead"
    ]), angular.module("ui.bootstrap.tpls", ["uib/template/accordion/accordion-group.html", "uib/template/accordion/accordion.html",
      "uib/template/alert/alert.html", "uib/template/carousel/carousel.html", "uib/template/carousel/slide.html",
      "uib/template/datepicker/datepicker.html", "uib/template/datepicker/day.html", "uib/template/datepicker/month.html",
      "uib/template/datepicker/year.html", "uib/template/datepickerPopup/popup.html", "uib/template/modal/window.html",
      "uib/template/pager/pager.html", "uib/template/pagination/pagination.html", "uib/template/tooltip/tooltip-html-popup.html",
      "uib/template/tooltip/tooltip-popup.html", "uib/template/tooltip/tooltip-template-popup.html", "uib/template/popover/popover-html.html",
      "uib/template/popover/popover-template.html", "uib/template/popover/popover.html", "uib/template/progressbar/bar.html",
      "uib/template/progressbar/progress.html", "uib/template/progressbar/progressbar.html", "uib/template/rating/rating.html",
      "uib/template/tabs/tab.html", "uib/template/tabs/tabset.html", "uib/template/timepicker/timepicker.html",
      "uib/template/typeahead/typeahead-match.html", "uib/template/typeahead/typeahead-popup.html"
    ]), angular.module("ui.bootstrap.collapse", []).directive("uibCollapse", ["$animate", "$q", "$parse", "$injector", function(e, t, n, r) {
      var i = r.has("$animateCss") ? r.get("$animateCss") : null;
      return {
        link: function(r, o, a) {
          var s = n(a.expanding),
            u = n(a.expanded),
            l = n(a.collapsing),
            c = n(a.collapsed),
            d = !1,
            p = {},
            f = {};

          function h(e) {
            return d ? {
              width: e.scrollWidth + "px"
            } : {
              height: e.scrollHeight + "px"
            }
          }

          function g() {
            o.removeClass("collapsing").addClass("collapse").css(p), u(r)
          }

          function m() {
            o.css(f), o.removeClass("collapsing").addClass("collapse"), c(r)
          }! function() {
            (d = !!("horizontal" in a)) ? (p = {
              width: ""
            }, f = {
              width: "0"
            }) : (p = {
              height: ""
            }, f = {
              height: "0"
            });
            r.$eval(a.uibCollapse) || o.addClass("in").addClass("collapse").attr("aria-expanded", !0).attr("aria-hidden", !1).css(p)
          }(), r.$watch(a.uibCollapse, (function(n) {
            n ? function() {
              if (!o.hasClass("collapse") && !o.hasClass("in")) return m();
              t.resolve(l(r)).then((function() {
                o.css(h(o[0])).removeClass("collapse").addClass("collapsing").attr("aria-expanded", !1).attr("aria-hidden", !0), i ?
                  i(o, {
                    removeClass: "in",
                    to: f
                  }).start().finally(m) : e.removeClass(o, "in", {
                    to: f
                  }).then(m)
              }), angular.noop)
            }() : o.hasClass("collapse") && o.hasClass("in") || t.resolve(s(r)).then((function() {
              o.removeClass("collapse").addClass("collapsing").attr("aria-expanded", !0).attr("aria-hidden", !1), i ? i(o, {
                addClass: "in",
                easing: "ease",
                css: {
                  overflow: "hidden"
                },
                to: h(o[0])
              }).start().finally(g) : e.addClass(o, "in", {
                css: {
                  overflow: "hidden"
                },
                to: h(o[0])
              }).then(g)
            }), angular.noop)
          }))
        }
      }
    }]), angular.module("ui.bootstrap.tabindex", []).directive("uibTabindexToggle", (function() {
      return {
        restrict: "A",
        link: function(e, t, n) {
          n.$observe("disabled", (function(e) {
            n.$set("tabindex", e ? -1 : null)
          }))
        }
      }
    })), angular.module("ui.bootstrap.accordion", ["ui.bootstrap.collapse", "ui.bootstrap.tabindex"]).constant("uibAccordionConfig", {
      closeOthers: !0
    }).controller("UibAccordionController", ["$scope", "$attrs", "uibAccordionConfig", function(e, t, n) {
      this.groups = [], this.closeOthers = function(r) {
        (angular.isDefined(t.closeOthers) ? e.$eval(t.closeOthers) : n.closeOthers) && angular.forEach(this.groups, (function(e) {
          e !== r && (e.isOpen = !1)
        }))
      }, this.addGroup = function(e) {
        var t = this;
        this.groups.push(e), e.$on("$destroy", (function(n) {
          t.removeGroup(e)
        }))
      }, this.removeGroup = function(e) {
        var t = this.groups.indexOf(e); - 1 !== t && this.groups.splice(t, 1)
      }
    }]).directive("uibAccordion", (function() {
      return {
        controller: "UibAccordionController",
        controllerAs: "accordion",
        transclude: !0,
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/accordion/accordion.html"
        }
      }
    })).directive("uibAccordionGroup", (function() {
      return {
        require: "^uibAccordion",
        transclude: !0,
        restrict: "A",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/accordion/accordion-group.html"
        },
        scope: {
          heading: "@",
          panelClass: "@?",
          isOpen: "=?",
          isDisabled: "=?"
        },
        controller: function() {
          this.setHeading = function(e) {
            this.heading = e
          }
        },
        link: function(e, t, n, r) {
          t.addClass("panel"), r.addGroup(e), e.openClass = n.openClass || "panel-open", e.panelClass = n.panelClass || "panel-default", e.$watch(
            "isOpen", (function(n) {
              t.toggleClass(e.openClass, !!n), n && r.closeOthers(e)
            })), e.toggleOpen = function(t) {
            e.isDisabled || t && 32 !== t.which || (e.isOpen = !e.isOpen)
          };
          var i = "accordiongroup-" + e.$id + "-" + Math.floor(1e4 * Math.random());
          e.headingId = i + "-tab", e.panelId = i + "-panel"
        }
      }
    })).directive("uibAccordionHeading", (function() {
      return {
        transclude: !0,
        template: "",
        replace: !0,
        require: "^uibAccordionGroup",
        link: function(e, t, n, r, i) {
          r.setHeading(i(e, angular.noop))
        }
      }
    })).directive("uibAccordionTransclude", (function() {
      return {
        require: "^uibAccordionGroup",
        link: function(e, t, n, r) {
          e.$watch((function() {
            return r[n.uibAccordionTransclude]
          }), (function(e) {
            if (e) {
              var n = angular.element(t[0].querySelector(
                "uib-accordion-header,data-uib-accordion-header,x-uib-accordion-header,uib\\:accordion-header,[uib-accordion-header],[data-uib-accordion-header],[x-uib-accordion-header]"
                ));
              n.html(""), n.append(e)
            }
          }))
        }
      }
    })), angular.module("ui.bootstrap.alert", []).controller("UibAlertController", ["$scope", "$element", "$attrs", "$interpolate", "$timeout",
      function(e, t, n, r, i) {
        e.closeable = !!n.close, t.addClass("alert"), n.$set("role", "alert"), e.closeable && t.addClass("alert-dismissible");
        var o = angular.isDefined(n.dismissOnTimeout) ? r(n.dismissOnTimeout)(e.$parent) : null;
        o && i((function() {
          e.close()
        }), parseInt(o, 10))
      }
    ]).directive("uibAlert", (function() {
      return {
        controller: "UibAlertController",
        controllerAs: "alert",
        restrict: "A",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/alert/alert.html"
        },
        transclude: !0,
        scope: {
          close: "&"
        }
      }
    })), angular.module("ui.bootstrap.buttons", []).constant("uibButtonConfig", {
      activeClass: "active",
      toggleEvent: "click"
    }).controller("UibButtonsController", ["uibButtonConfig", function(e) {
      this.activeClass = e.activeClass || "active", this.toggleEvent = e.toggleEvent || "click"
    }]).directive("uibBtnRadio", ["$parse", function(e) {
      return {
        require: ["uibBtnRadio", "ngModel"],
        controller: "UibButtonsController",
        controllerAs: "buttons",
        link: function(t, n, r, i) {
          var o = i[0],
            a = i[1],
            s = e(r.uibUncheckable);
          n.find("input").css({
            display: "none"
          }), a.$render = function() {
            n.toggleClass(o.activeClass, angular.equals(a.$modelValue, t.$eval(r.uibBtnRadio)))
          }, n.on(o.toggleEvent, (function() {
            if (!r.disabled) {
              var e = n.hasClass(o.activeClass);
              e && !angular.isDefined(r.uncheckable) || t.$apply((function() {
                a.$setViewValue(e ? null : t.$eval(r.uibBtnRadio)), a.$render()
              }))
            }
          })), r.uibUncheckable && t.$watch(s, (function(e) {
            r.$set("uncheckable", e ? "" : void 0)
          }))
        }
      }
    }]).directive("uibBtnCheckbox", (function() {
      return {
        require: ["uibBtnCheckbox", "ngModel"],
        controller: "UibButtonsController",
        controllerAs: "button",
        link: function(e, t, n, r) {
          var i = r[0],
            o = r[1];

          function a() {
            return s(n.btnCheckboxTrue, !0)
          }

          function s(t, n) {
            return angular.isDefined(t) ? e.$eval(t) : n
          }
          t.find("input").css({
            display: "none"
          }), o.$render = function() {
            t.toggleClass(i.activeClass, angular.equals(o.$modelValue, a()))
          }, t.on(i.toggleEvent, (function() {
            n.disabled || e.$apply((function() {
              o.$setViewValue(t.hasClass(i.activeClass) ? s(n.btnCheckboxFalse, !1) : a()), o.$render()
            }))
          }))
        }
      }
    })), angular.module("ui.bootstrap.carousel", []).controller("UibCarouselController", ["$scope", "$element", "$interval", "$timeout", "$animate",
      function(e, t, n, r, i) {
        var o, a, s = this,
          u = s.slides = e.slides = [],
          l = "uib-slideDirection",
          c = e.active,
          d = !1;

        function p(e) {
          for (var t = 0; t < u.length; t++) u[t].slide.active = t === e
        }

        function f(e) {
          for (var t = 0; t < u.length; t++)
            if (u[t].slide === e) return t
        }

        function h() {
          o && (n.cancel(o), o = null)
        }

        function g() {
          h();
          var t = +e.interval;
          !isNaN(t) && t > 0 && (o = n(m, t))
        }

        function m() {
          var t = +e.interval;
          a && !isNaN(t) && t > 0 && u.length ? e.next() : e.pause()
        }
        t.addClass("carousel"), s.addSlide = function(t, n) {
          u.push({
            slide: t,
            element: n
          }), u.sort((function(e, t) {
            return +e.slide.index - +t.slide.index
          })), (t.index === e.active || 1 === u.length && !angular.isNumber(e.active)) && (e.$currentTransition && (e.$currentTransition = null),
            c = t.index, e.active = t.index, p(c), s.select(u[f(t)]), 1 === u.length && e.play())
        }, s.getCurrentIndex = function() {
          for (var e = 0; e < u.length; e++)
            if (u[e].slide.index === c) return e
        }, s.next = e.next = function() {
          var t = (s.getCurrentIndex() + 1) % u.length;
          if (0 !== t || !e.noWrap()) return s.select(u[t], "next");
          e.pause()
        }, s.prev = e.prev = function() {
          var t = s.getCurrentIndex() - 1 < 0 ? u.length - 1 : s.getCurrentIndex() - 1;
          if (!e.noWrap() || t !== u.length - 1) return s.select(u[t], "prev");
          e.pause()
        }, s.removeSlide = function(t) {
          var n = f(t);
          u.splice(n, 1), u.length > 0 && c === n ? n >= u.length ? (c = u.length - 1, e.active = c, p(c), s.select(u[u.length - 1])) : (c = n, e
            .active = c, p(c), s.select(u[n])) : c > n && (c--, e.active = c), 0 === u.length && (c = null, e.active = null)
        }, s.select = e.select = function(n, r) {
          var o = f(n.slide);
          void 0 === r && (r = o > s.getCurrentIndex() ? "next" : "prev"), n.slide.index === c || e.$currentTransition || function(n, r, o) {
            if (d) return;
            if (angular.extend(n, {
                direction: o
              }), angular.extend(u[c].slide || {}, {
                direction: o
              }), i.enabled(t) && !e.$currentTransition && u[r].element && s.slides.length > 1) {
              u[r].element.data(l, n.direction);
              var a = s.getCurrentIndex();
              angular.isNumber(a) && u[a].element && u[a].element.data(l, n.direction), e.$currentTransition = !0, i.on("addClass", u[r].element,
                (function(t, n) {
                  "close" === n && (e.$currentTransition = null, i.off("addClass", t))
                }))
            }
            e.active = n.index, c = n.index, p(r), g()
          }(n.slide, o, r)
        }, e.indexOfSlide = function(e) {
          return +e.slide.index
        }, e.isActive = function(t) {
          return e.active === t.slide.index
        }, e.isPrevDisabled = function() {
          return 0 === e.active && e.noWrap()
        }, e.isNextDisabled = function() {
          return e.active === u.length - 1 && e.noWrap()
        }, e.pause = function() {
          e.noPause || (a = !1, h())
        }, e.play = function() {
          a || (a = !0, g())
        }, t.on("mouseenter", e.pause), t.on("mouseleave", e.play), e.$on("$destroy", (function() {
          d = !0, h()
        })), e.$watch("noTransition", (function(e) {
          i.enabled(t, !e)
        })), e.$watch("interval", g), e.$watchCollection("slides", (function(t) {
          t.length || (e.$currentTransition = null)
        })), e.$watch("active", (function(e) {
          if (angular.isNumber(e) && c !== e) {
            for (var t = 0; t < u.length; t++)
              if (u[t].slide.index === e) {
                e = t;
                break
              } u[e] && (p(e), s.select(u[e]), c = e)
          }
        }))
      }
    ]).directive("uibCarousel", (function() {
      return {
        transclude: !0,
        controller: "UibCarouselController",
        controllerAs: "carousel",
        restrict: "A",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/carousel/carousel.html"
        },
        scope: {
          active: "=",
          interval: "=",
          noTransition: "=",
          noPause: "=",
          noWrap: "&"
        }
      }
    })).directive("uibSlide", ["$animate", function(e) {
      return {
        require: "^uibCarousel",
        restrict: "A",
        transclude: !0,
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/carousel/slide.html"
        },
        scope: {
          actual: "=?",
          index: "=?"
        },
        link: function(t, n, r, i) {
          n.addClass("item"), i.addSlide(t, n), t.$on("$destroy", (function() {
            i.removeSlide(t)
          })), t.$watch("active", (function(t) {
            e[t ? "addClass" : "removeClass"](n, "active")
          }))
        }
      }
    }]).animation(".item", ["$animateCss", function(e) {
      var t = "uib-slideDirection";

      function n(e, t, n) {
        e.removeClass(t), n && n()
      }
      return {
        beforeAddClass: function(r, i, o) {
          if ("active" === i) {
            var a = r.data(t),
              s = "next" === a ? "left" : "right",
              u = n.bind(this, r, s + " " + a, o);
            return r.addClass(a), e(r, {
                addClass: s
              }).start().done(u),
              function() {
                !0
              }
          }
          o()
        },
        beforeRemoveClass: function(r, i, o) {
          if ("active" === i) {
            var a = "next" === r.data(t) ? "left" : "right",
              s = n.bind(this, r, a, o);
            return e(r, {
                addClass: a
              }).start().done(s),
              function() {
                !0
              }
          }
          o()
        }
      }
    }]), angular.module("ui.bootstrap.dateparser", []).service("uibDateParser", ["$log", "$locale", "dateFilter", "orderByFilter", "filterFilter",
      function(e, t, n, r, i) {
        var o, a, s = /[\\\^\$\*\+\?\|\[\]\(\)\.\{\}]/g;

        function u(e) {
          return i(a, {
            key: e
          }, !0)[0]
        }

        function l(e, t, n) {
          return function() {
            return e.substr(t + 1, n - t - 1)
          }
        }

        function c(e, t) {
          for (var n = e.substr(t), r = 0; r < a.length; r++)
            if (new RegExp("^" + a[r].key).test(n)) {
              var i = a[r];
              return {
                endIdx: t + i.key.length,
                parser: i.formatter
              }
            } return {
            endIdx: t + 1,
            parser: function() {
              return n.charAt(0)
            }
          }
        }

        function d(e) {
          return parseInt(e, 10)
        }

        function p(e, t) {
          e = e.replace(/:/g, "");
          var n = Date.parse("Jan 01, 1970 00:00:00 " + e) / 6e4;
          return isNaN(n) ? t : n
        }

        function f(e, t) {
          return (e = new Date(e.getTime())).setMinutes(e.getMinutes() + t), e
        }

        function h(e, t, n) {
          n = n ? -1 : 1;
          var r = e.getTimezoneOffset();
          return f(e, n * (p(t, r) - r))
        }
        this.init = function() {
          o = t.id, this.parsers = {}, this.formatters = {}, a = [{
            key: "yyyy",
            regex: "\\d{4}",
            apply: function(e) {
              this.year = +e
            },
            formatter: function(e) {
              var t = new Date;
              return t.setFullYear(Math.abs(e.getFullYear())), n(t, "yyyy")
            }
          }, {
            key: "yy",
            regex: "\\d{2}",
            apply: function(e) {
              e = +e, this.year = e < 69 ? e + 2e3 : e + 1900
            },
            formatter: function(e) {
              var t = new Date;
              return t.setFullYear(Math.abs(e.getFullYear())), n(t, "yy")
            }
          }, {
            key: "y",
            regex: "\\d{1,4}",
            apply: function(e) {
              this.year = +e
            },
            formatter: function(e) {
              var t = new Date;
              return t.setFullYear(Math.abs(e.getFullYear())), n(t, "y")
            }
          }, {
            key: "M!",
            regex: "0?[1-9]|1[0-2]",
            apply: function(e) {
              this.month = e - 1
            },
            formatter: function(e) {
              var t = e.getMonth();
              return /^[0-9]$/.test(t) ? n(e, "MM") : n(e, "M")
            }
          }, {
            key: "MMMM",
            regex: t.DATETIME_FORMATS.MONTH.join("|"),
            apply: function(e) {
              this.month = t.DATETIME_FORMATS.MONTH.indexOf(e)
            },
            formatter: function(e) {
              return n(e, "MMMM")
            }
          }, {
            key: "MMM",
            regex: t.DATETIME_FORMATS.SHORTMONTH.join("|"),
            apply: function(e) {
              this.month = t.DATETIME_FORMATS.SHORTMONTH.indexOf(e)
            },
            formatter: function(e) {
              return n(e, "MMM")
            }
          }, {
            key: "MM",
            regex: "0[1-9]|1[0-2]",
            apply: function(e) {
              this.month = e - 1
            },
            formatter: function(e) {
              return n(e, "MM")
            }
          }, {
            key: "M",
            regex: "[1-9]|1[0-2]",
            apply: function(e) {
              this.month = e - 1
            },
            formatter: function(e) {
              return n(e, "M")
            }
          }, {
            key: "d!",
            regex: "[0-2]?[0-9]{1}|3[0-1]{1}",
            apply: function(e) {
              this.date = +e
            },
            formatter: function(e) {
              var t = e.getDate();
              return /^[1-9]$/.test(t) ? n(e, "dd") : n(e, "d")
            }
          }, {
            key: "dd",
            regex: "[0-2][0-9]{1}|3[0-1]{1}",
            apply: function(e) {
              this.date = +e
            },
            formatter: function(e) {
              return n(e, "dd")
            }
          }, {
            key: "d",
            regex: "[1-2]?[0-9]{1}|3[0-1]{1}",
            apply: function(e) {
              this.date = +e
            },
            formatter: function(e) {
              return n(e, "d")
            }
          }, {
            key: "EEEE",
            regex: t.DATETIME_FORMATS.DAY.join("|"),
            formatter: function(e) {
              return n(e, "EEEE")
            }
          }, {
            key: "EEE",
            regex: t.DATETIME_FORMATS.SHORTDAY.join("|"),
            formatter: function(e) {
              return n(e, "EEE")
            }
          }, {
            key: "HH",
            regex: "(?:0|1)[0-9]|2[0-3]",
            apply: function(e) {
              this.hours = +e
            },
            formatter: function(e) {
              return n(e, "HH")
            }
          }, {
            key: "hh",
            regex: "0[0-9]|1[0-2]",
            apply: function(e) {
              this.hours = +e
            },
            formatter: function(e) {
              return n(e, "hh")
            }
          }, {
            key: "H",
            regex: "1?[0-9]|2[0-3]",
            apply: function(e) {
              this.hours = +e
            },
            formatter: function(e) {
              return n(e, "H")
            }
          }, {
            key: "h",
            regex: "[0-9]|1[0-2]",
            apply: function(e) {
              this.hours = +e
            },
            formatter: function(e) {
              return n(e, "h")
            }
          }, {
            key: "mm",
            regex: "[0-5][0-9]",
            apply: function(e) {
              this.minutes = +e
            },
            formatter: function(e) {
              return n(e, "mm")
            }
          }, {
            key: "m",
            regex: "[0-9]|[1-5][0-9]",
            apply: function(e) {
              this.minutes = +e
            },
            formatter: function(e) {
              return n(e, "m")
            }
          }, {
            key: "sss",
            regex: "[0-9][0-9][0-9]",
            apply: function(e) {
              this.milliseconds = +e
            },
            formatter: function(e) {
              return n(e, "sss")
            }
          }, {
            key: "ss",
            regex: "[0-5][0-9]",
            apply: function(e) {
              this.seconds = +e
            },
            formatter: function(e) {
              return n(e, "ss")
            }
          }, {
            key: "s",
            regex: "[0-9]|[1-5][0-9]",
            apply: function(e) {
              this.seconds = +e
            },
            formatter: function(e) {
              return n(e, "s")
            }
          }, {
            key: "a",
            regex: t.DATETIME_FORMATS.AMPMS.join("|"),
            apply: function(e) {
              12 === this.hours && (this.hours = 0), "PM" === e && (this.hours += 12)
            },
            formatter: function(e) {
              return n(e, "a")
            }
          }, {
            key: "Z",
            regex: "[+-]\\d{4}",
            apply: function(e) {
              var t = e.match(/([+-])(\d{2})(\d{2})/),
                n = t[1],
                r = t[2],
                i = t[3];
              this.hours += d(n + r), this.minutes += d(n + i)
            },
            formatter: function(e) {
              return n(e, "Z")
            }
          }, {
            key: "ww",
            regex: "[0-4][0-9]|5[0-3]",
            formatter: function(e) {
              return n(e, "ww")
            }
          }, {
            key: "w",
            regex: "[0-9]|[1-4][0-9]|5[0-3]",
            formatter: function(e) {
              return n(e, "w")
            }
          }, {
            key: "GGGG",
            regex: t.DATETIME_FORMATS.ERANAMES.join("|").replace(/\s/g, "\\s"),
            formatter: function(e) {
              return n(e, "GGGG")
            }
          }, {
            key: "GGG",
            regex: t.DATETIME_FORMATS.ERAS.join("|"),
            formatter: function(e) {
              return n(e, "GGG")
            }
          }, {
            key: "GG",
            regex: t.DATETIME_FORMATS.ERAS.join("|"),
            formatter: function(e) {
              return n(e, "GG")
            }
          }, {
            key: "G",
            regex: t.DATETIME_FORMATS.ERAS.join("|"),
            formatter: function(e) {
              return n(e, "G")
            }
          }], angular.version.major >= 1 && angular.version.minor > 4 && a.push({
            key: "LLLL",
            regex: t.DATETIME_FORMATS.STANDALONEMONTH.join("|"),
            apply: function(e) {
              this.month = t.DATETIME_FORMATS.STANDALONEMONTH.indexOf(e)
            },
            formatter: function(e) {
              return n(e, "LLLL")
            }
          })
        }, this.init(), this.getParser = function(e) {
          var t = u(e);
          return t && t.apply || null
        }, this.overrideParser = function(e, t) {
          var n = u(e);
          n && angular.isFunction(t) && (this.parsers = {}, n.apply = t)
        }.bind(this), this.filter = function(e, n) {
          return angular.isDate(e) && !isNaN(e) && n ? (n = t.DATETIME_FORMATS[n] || n, t.id !== o && this.init(), this.formatters[n] || (this
            .formatters[n] = function(e) {
              for (var t, n, r = [], i = 0; i < e.length;)
                if (angular.isNumber(n)) {
                  if ("'" === e.charAt(i))(i + 1 >= e.length || "'" !== e.charAt(i + 1)) && (r.push(l(e, n, i)), n = null);
                  else if (i === e.length)
                    for (; n < e.length;) t = c(e, n), r.push(t), n = t.endIdx;
                  i++
                } else "'" !== e.charAt(i) ? (t = c(e, i), r.push(t.parser), i = t.endIdx) : (n = i, i++);
              return r
            }(n)), this.formatters[n].reduce((function(t, n) {
            return t + n(e)
          }), "")) : ""
        }, this.parse = function(n, i, u) {
          if (!angular.isString(n) || !i) return n;
          i = (i = t.DATETIME_FORMATS[i] || i).replace(s, "\\$&"), t.id !== o && this.init(), this.parsers[i] || (this.parsers[i] = function(e) {
            var t = [],
              n = e.split(""),
              i = e.indexOf("'");
            if (i > -1) {
              var o = !1;
              e = e.split("");
              for (var s = i; s < e.length; s++) o ? ("'" === e[s] && (s + 1 < e.length && "'" === e[s + 1] ? (e[s + 1] = "$", n[s + 1] = "") :
                (n[s] = "", o = !1)), e[s] = "$") : "'" === e[s] && (e[s] = "$", n[s] = "", o = !0);
              e = e.join("")
            }
            return angular.forEach(a, (function(r) {
              var i = e.indexOf(r.key);
              if (i > -1) {
                e = e.split(""), n[i] = "(" + r.regex + ")", e[i] = "$";
                for (var o = i + 1, a = i + r.key.length; o < a; o++) n[o] = "", e[o] = "$";
                e = e.join(""), t.push({
                  index: i,
                  key: r.key,
                  apply: r.apply,
                  matcher: r.regex
                })
              }
            })), {
              regex: new RegExp("^" + n.join("") + "$"),
              map: r(t, "index")
            }
          }(i));
          var l = this.parsers[i],
            c = l.regex,
            d = l.map,
            p = n.match(c),
            f = !1;
          if (p && p.length) {
            var h, g;
            angular.isDate(u) && !isNaN(u.getTime()) ? h = {
              year: u.getFullYear(),
              month: u.getMonth(),
              date: u.getDate(),
              hours: u.getHours(),
              minutes: u.getMinutes(),
              seconds: u.getSeconds(),
              milliseconds: u.getMilliseconds()
            } : (u && e.warn("dateparser:", "baseDate is not a valid date"), h = {
              year: 1900,
              month: 0,
              date: 1,
              hours: 0,
              minutes: 0,
              seconds: 0,
              milliseconds: 0
            });
            for (var m = 1, v = p.length; m < v; m++) {
              var $ = d[m - 1];
              "Z" === $.matcher && (f = !0), $.apply && $.apply.call(h, p[m])
            }
            var b = f ? Date.prototype.setUTCFullYear : Date.prototype.setFullYear,
              y = f ? Date.prototype.setUTCHours : Date.prototype.setHours;
            return function(e, t, n) {
              if (n < 1) return !1;
              if (1 === t && n > 28) return 29 === n && (e % 4 == 0 && e % 100 != 0 || e % 400 == 0);
              if (3 === t || 5 === t || 8 === t || 10 === t) return n < 31;
              return !0
            }(h.year, h.month, h.date) && (!angular.isDate(u) || isNaN(u.getTime()) || f ? (g = new Date(0), b.call(g, h.year, h.month, h.date), y
              .call(g, h.hours || 0, h.minutes || 0, h.seconds || 0, h.milliseconds || 0)) : (g = new Date(u), b.call(g, h.year, h.month, h
              .date), y.call(g, h.hours, h.minutes, h.seconds, h.milliseconds))), g
          }
        }, this.toTimezone = function(e, t) {
          return e && t ? h(e, t) : e
        }, this.fromTimezone = function(e, t) {
          return e && t ? h(e, t, !0) : e
        }, this.timezoneToOffset = p, this.addDateMinutes = f, this.convertTimezoneToLocal = h
      }
    ]), angular.module("ui.bootstrap.isClass", []).directive("uibIsClass", ["$animate", function(e) {
      var t = /^\s*([\s\S]+?)\s+on\s+([\s\S]+?)\s*$/,
        n = /^\s*([\s\S]+?)\s+for\s+([\s\S]+?)\s*$/;
      return {
        restrict: "A",
        compile: function(r, i) {
          var o = [],
            a = [],
            s = {},
            u = i.uibIsClass.match(t),
            l = u[2],
            c = u[1].split(",");
          return function(t, r, i) {
            o.push(t), a.push({
              scope: t,
              element: r
            }), c.forEach((function(r, i) {
              ! function(t, r) {
                var i = t.match(n),
                  o = r.$eval(i[1]),
                  u = i[2],
                  c = s[t];
                if (!c) {
                  var d = function(t) {
                    var n = null;
                    a.some((function(e) {
                      if (e.scope.$eval(l) === t) return n = e, !0
                    })), c.lastActivated !== n && (c.lastActivated && e.removeClass(c.lastActivated.element, o), n && e.addClass(n
                      .element, o), c.lastActivated = n)
                  };
                  s[t] = c = {
                    lastActivated: null,
                    scope: r,
                    watchFn: d,
                    compareWithExp: u,
                    watcher: r.$watch(u, d)
                  }
                }
                c.watchFn(r.$eval(u))
              }(r, t)
            })), t.$on("$destroy", d)
          };

          function d(e) {
            var t = e.targetScope,
              n = o.indexOf(t);
            if (o.splice(n, 1), a.splice(n, 1), o.length) {
              var r = o[0];
              angular.forEach(s, (function(e) {
                e.scope === t && (e.watcher = r.$watch(e.compareWithExp, e.watchFn), e.scope = r)
              }))
            } else s = {}
          }
        }
      }
    }]), angular.module("ui.bootstrap.datepicker", ["ui.bootstrap.dateparser", "ui.bootstrap.isClass"]).value("$datepickerSuppressError", !1).value(
      "$datepickerLiteralWarning", !0).constant("uibDatepickerConfig", {
      datepickerMode: "day",
      formatDay: "dd",
      formatMonth: "MMMM",
      formatYear: "yyyy",
      formatDayHeader: "EEE",
      formatDayTitle: "MMMM yyyy",
      formatMonthTitle: "yyyy",
      maxDate: null,
      maxMode: "year",
      minDate: null,
      minMode: "day",
      monthColumns: 3,
      ngModelOptions: {},
      shortcutPropagation: !1,
      showWeeks: !0,
      yearColumns: 5,
      yearRows: 4
    }).controller("UibDatepickerController", ["$scope", "$element", "$attrs", "$parse", "$interpolate", "$locale", "$log", "dateFilter",
      "uibDatepickerConfig", "$datepickerLiteralWarning", "$datepickerSuppressError", "uibDateParser",
      function(e, t, n, r, i, o, a, s, u, l, c, d) {
        var p = this,
          f = {
            $setViewValue: angular.noop
          },
          h = {},
          g = [];
        t.addClass("uib-datepicker"), n.$set("role", "application"), e.datepickerOptions || (e.datepickerOptions = {}), this.modes = ["day", "month",
            "year"
          ], ["customClass", "dateDisabled", "datepickerMode", "formatDay", "formatDayHeader", "formatDayTitle", "formatMonth", "formatMonthTitle",
            "formatYear", "maxDate", "maxMode", "minDate", "minMode", "monthColumns", "showWeeks", "shortcutPropagation", "startingDay",
            "yearColumns", "yearRows"
          ].forEach((function(t) {
            switch (t) {
              case "customClass":
              case "dateDisabled":
                e[t] = e.datepickerOptions[t] || angular.noop;
                break;
              case "datepickerMode":
                e.datepickerMode = angular.isDefined(e.datepickerOptions.datepickerMode) ? e.datepickerOptions.datepickerMode : u.datepickerMode;
                break;
              case "formatDay":
              case "formatDayHeader":
              case "formatDayTitle":
              case "formatMonth":
              case "formatMonthTitle":
              case "formatYear":
                p[t] = angular.isDefined(e.datepickerOptions[t]) ? i(e.datepickerOptions[t])(e.$parent) : u[t];
                break;
              case "monthColumns":
              case "showWeeks":
              case "shortcutPropagation":
              case "yearColumns":
              case "yearRows":
                p[t] = angular.isDefined(e.datepickerOptions[t]) ? e.datepickerOptions[t] : u[t];
                break;
              case "startingDay":
                angular.isDefined(e.datepickerOptions.startingDay) ? p.startingDay = e.datepickerOptions.startingDay : angular.isNumber(u
                  .startingDay) ? p.startingDay = u.startingDay : p.startingDay = (o.DATETIME_FORMATS.FIRSTDAYOFWEEK + 8) % 7;
                break;
              case "maxDate":
              case "minDate":
                e.$watch("datepickerOptions." + t, (function(e) {
                  e ? angular.isDate(e) ? p[t] = d.fromTimezone(new Date(e), h.getOption("timezone")) : (l && a.warn(
                    "Literal date support has been deprecated, please switch to date object usage"), p[t] = new Date(s(e, "medium"))) : p[
                    t] = u[t] ? d.fromTimezone(new Date(u[t]), h.getOption("timezone")) : null, p.refreshView()
                }));
                break;
              case "maxMode":
              case "minMode":
                e.datepickerOptions[t] ? e.$watch((function() {
                  return e.datepickerOptions[t]
                }), (function(n) {
                  p[t] = e[t] = angular.isDefined(n) ? n : e.datepickerOptions[t], ("minMode" === t && p.modes.indexOf(e.datepickerOptions
                      .datepickerMode) < p.modes.indexOf(p[t]) || "maxMode" === t && p.modes.indexOf(e.datepickerOptions.datepickerMode) >
                    p.modes.indexOf(p[t])) && (e.datepickerMode = p[t], e.datepickerOptions.datepickerMode = p[t])
                })) : p[t] = e[t] = u[t] || null
            }
          })), e.uniqueId = "datepicker-" + e.$id + "-" + Math.floor(1e4 * Math.random()), e.disabled = angular.isDefined(n.disabled) || !1, angular
          .isDefined(n.ngDisabled) && g.push(e.$parent.$watch(n.ngDisabled, (function(t) {
            e.disabled = t, p.refreshView()
          }))), e.isActive = function(t) {
            return 0 === p.compare(t.date, p.activeDate) && (e.activeDateId = t.uid, !0)
          }, this.init = function(t) {
            h = function(t) {
              var n;
              if (angular.version.minor < 6)(n = t.$options || e.datepickerOptions.ngModelOptions || u.ngModelOptions || {}).getOption = function(
              e) {
                return n[e]
              };
              else {
                var r = t.$options.getOption("timezone") || (e.datepickerOptions.ngModelOptions ? e.datepickerOptions.ngModelOptions.timezone :
                  null) || (u.ngModelOptions ? u.ngModelOptions.timezone : null);
                n = t.$options.createChild(u.ngModelOptions).createChild(e.datepickerOptions.ngModelOptions).createChild(t.$options).createChild({
                  timezone: r
                })
              }
              return n
            }(f = t), e.datepickerOptions.initDate ? (p.activeDate = d.fromTimezone(e.datepickerOptions.initDate, h.getOption("timezone")) ||
              new Date, e.$watch("datepickerOptions.initDate", (function(e) {
                e && (f.$isEmpty(f.$modelValue) || f.$invalid) && (p.activeDate = d.fromTimezone(e, h.getOption("timezone")), p.refreshView())
              }))) : p.activeDate = new Date;
            var n = f.$modelValue ? new Date(f.$modelValue) : new Date;
            this.activeDate = isNaN(n) ? d.fromTimezone(new Date, h.getOption("timezone")) : d.fromTimezone(n, h.getOption("timezone")), f.$render =
              function() {
                p.render()
              }
          }, this.render = function() {
            if (f.$viewValue) {
              var e = new Date(f.$viewValue);
              !isNaN(e) ? this.activeDate = d.fromTimezone(e, h.getOption("timezone")) : c || a.error(
                'Datepicker directive: "ng-model" value must be a Date object')
            }
            this.refreshView()
          }, this.refreshView = function() {
            if (this.element) {
              e.selectedDt = null, this._refreshView(), e.activeDt && (e.activeDateId = e.activeDt.uid);
              var t = f.$viewValue ? new Date(f.$viewValue) : null;
              t = d.fromTimezone(t, h.getOption("timezone")), f.$setValidity("dateDisabled", !t || this.element && !this.isDisabled(t))
            }
          }, this.createDateObject = function(t, n) {
            var r = f.$viewValue ? new Date(f.$viewValue) : null;
            r = d.fromTimezone(r, h.getOption("timezone"));
            var i = new Date;
            i = d.fromTimezone(i, h.getOption("timezone"));
            var o = this.compare(t, i),
              a = {
                date: t,
                label: d.filter(t, n),
                selected: r && 0 === this.compare(t, r),
                disabled: this.isDisabled(t),
                past: o < 0,
                current: 0 === o,
                future: o > 0,
                customClass: this.customClass(t) || null
              };
            return r && 0 === this.compare(t, r) && (e.selectedDt = a), p.activeDate && 0 === this.compare(a.date, p.activeDate) && (e.activeDt = a),
              a
          }, this.isDisabled = function(t) {
            return e.disabled || this.minDate && this.compare(t, this.minDate) < 0 || this.maxDate && this.compare(t, this.maxDate) > 0 || e
              .dateDisabled && e.dateDisabled({
                date: t,
                mode: e.datepickerMode
              })
          }, this.customClass = function(t) {
            return e.customClass({
              date: t,
              mode: e.datepickerMode
            })
          }, this.split = function(e, t) {
            for (var n = []; e.length > 0;) n.push(e.splice(0, t));
            return n
          }, e.select = function(t) {
            if (e.datepickerMode === p.minMode) {
              var n = f.$viewValue ? d.fromTimezone(new Date(f.$viewValue), h.getOption("timezone")) : new Date(0, 0, 0, 0, 0, 0, 0);
              n.setFullYear(t.getFullYear(), t.getMonth(), t.getDate()), n = d.toTimezone(n, h.getOption("timezone")), f.$setViewValue(n), f.$render()
            } else p.activeDate = t, m(p.modes[p.modes.indexOf(e.datepickerMode) - 1]), e.$emit("uib:datepicker.mode");
            e.$broadcast("uib:datepicker.focus")
          }, e.move = function(e) {
            var t = p.activeDate.getFullYear() + e * (p.step.years || 0),
              n = p.activeDate.getMonth() + e * (p.step.months || 0);
            p.activeDate.setFullYear(t, n, 1), p.refreshView()
          }, e.toggleMode = function(t) {
            t = t || 1, e.datepickerMode === p.maxMode && 1 === t || e.datepickerMode === p.minMode && -1 === t || (m(p.modes[p.modes.indexOf(e
              .datepickerMode) + t]), e.$emit("uib:datepicker.mode"))
          }, e.keys = {
            13: "enter",
            32: "space",
            33: "pageup",
            34: "pagedown",
            35: "end",
            36: "home",
            37: "left",
            38: "up",
            39: "right",
            40: "down"
          };

        function m(t) {
          e.datepickerMode = t, e.datepickerOptions.datepickerMode = t
        }
        e.$on("uib:datepicker.focus", (function() {
          p.element[0].focus()
        })), e.keydown = function(t) {
          var n = e.keys[t.which];
          if (n && !t.shiftKey && !t.altKey && !e.disabled)
            if (t.preventDefault(), p.shortcutPropagation || t.stopPropagation(), "enter" === n || "space" === n) {
              if (p.isDisabled(p.activeDate)) return;
              e.select(p.activeDate)
            } else !t.ctrlKey || "up" !== n && "down" !== n ? (p.handleKeyDown(n, t), p.refreshView()) : e.toggleMode("up" === n ? 1 : -1)
        }, t.on("keydown", (function(t) {
          e.$apply((function() {
            e.keydown(t)
          }))
        })), e.$on("$destroy", (function() {
          for (; g.length;) g.shift()()
        }))
      }
    ]).controller("UibDaypickerController", ["$scope", "$element", "dateFilter", function(e, t, n) {
      var r = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

      function i(e, t) {
        return 1 !== t || e % 4 != 0 || e % 100 == 0 && e % 400 != 0 ? r[t] : 29
      }

      function o(e) {
        var t = new Date(e);
        t.setDate(t.getDate() + 4 - (t.getDay() || 7));
        var n = t.getTime();
        return t.setMonth(0), t.setDate(1), Math.floor(Math.round((n - t) / 864e5) / 7) + 1
      }
      this.step = {
        months: 1
      }, this.element = t, this.init = function(t) {
        angular.extend(t, this), e.showWeeks = t.showWeeks, t.refreshView()
      }, this.getDates = function(e, t) {
        for (var n, r = new Array(t), i = new Date(e), o = 0; o < t;) n = new Date(i), r[o++] = n, i.setDate(i.getDate() + 1);
        return r
      }, this._refreshView = function() {
        var t = this.activeDate.getFullYear(),
          r = this.activeDate.getMonth(),
          i = new Date(this.activeDate);
        i.setFullYear(t, r, 1);
        var a = this.startingDay - i.getDay(),
          s = a > 0 ? 7 - a : -a,
          u = new Date(i);
        s > 0 && u.setDate(1 - s);
        for (var l = this.getDates(u, 42), c = 0; c < 42; c++) l[c] = angular.extend(this.createDateObject(l[c], this.formatDay), {
          secondary: l[c].getMonth() !== r,
          uid: e.uniqueId + "-" + c
        });
        e.labels = new Array(7);
        for (var d = 0; d < 7; d++) e.labels[d] = {
          abbr: n(l[d].date, this.formatDayHeader),
          full: n(l[d].date, "EEEE")
        };
        if (e.title = n(this.activeDate, this.formatDayTitle), e.rows = this.split(l, 7), e.showWeeks) {
          e.weekNumbers = [];
          for (var p = (11 - this.startingDay) % 7, f = e.rows.length, h = 0; h < f; h++) e.weekNumbers.push(o(e.rows[h][p].date))
        }
      }, this.compare = function(e, t) {
        var n = new Date(e.getFullYear(), e.getMonth(), e.getDate()),
          r = new Date(t.getFullYear(), t.getMonth(), t.getDate());
        return n.setFullYear(e.getFullYear()), r.setFullYear(t.getFullYear()), n - r
      }, this.handleKeyDown = function(e, t) {
        var n = this.activeDate.getDate();
        if ("left" === e) n -= 1;
        else if ("up" === e) n -= 7;
        else if ("right" === e) n += 1;
        else if ("down" === e) n += 7;
        else if ("pageup" === e || "pagedown" === e) {
          var r = this.activeDate.getMonth() + ("pageup" === e ? -1 : 1);
          this.activeDate.setMonth(r, 1), n = Math.min(i(this.activeDate.getFullYear(), this.activeDate.getMonth()), n)
        } else "home" === e ? n = 1 : "end" === e && (n = i(this.activeDate.getFullYear(), this.activeDate.getMonth()));
        this.activeDate.setDate(n)
      }
    }]).controller("UibMonthpickerController", ["$scope", "$element", "dateFilter", function(e, t, n) {
      this.step = {
        years: 1
      }, this.element = t, this.init = function(e) {
        angular.extend(e, this), e.refreshView()
      }, this._refreshView = function() {
        for (var t, r = new Array(12), i = this.activeDate.getFullYear(), o = 0; o < 12; o++)(t = new Date(this.activeDate)).setFullYear(i, o, 1),
          r[o] = angular.extend(this.createDateObject(t, this.formatMonth), {
            uid: e.uniqueId + "-" + o
          });
        e.title = n(this.activeDate, this.formatMonthTitle), e.rows = this.split(r, this.monthColumns), e.yearHeaderColspan = this.monthColumns >
          3 ? this.monthColumns - 2 : 1
      }, this.compare = function(e, t) {
        var n = new Date(e.getFullYear(), e.getMonth()),
          r = new Date(t.getFullYear(), t.getMonth());
        return n.setFullYear(e.getFullYear()), r.setFullYear(t.getFullYear()), n - r
      }, this.handleKeyDown = function(e, t) {
        var n = this.activeDate.getMonth();
        if ("left" === e) n -= 1;
        else if ("up" === e) n -= this.monthColumns;
        else if ("right" === e) n += 1;
        else if ("down" === e) n += this.monthColumns;
        else if ("pageup" === e || "pagedown" === e) {
          var r = this.activeDate.getFullYear() + ("pageup" === e ? -1 : 1);
          this.activeDate.setFullYear(r)
        } else "home" === e ? n = 0 : "end" === e && (n = 11);
        this.activeDate.setMonth(n)
      }
    }]).controller("UibYearpickerController", ["$scope", "$element", "dateFilter", function(e, t, n) {
      var r, i;

      function o(e) {
        return parseInt((e - 1) / i, 10) * i + 1
      }
      this.element = t, this.yearpickerInit = function() {
        r = this.yearColumns, i = this.yearRows * r, this.step = {
          years: i
        }
      }, this._refreshView = function() {
        for (var t, n = new Array(i), a = 0, s = o(this.activeDate.getFullYear()); a < i; a++)(t = new Date(this.activeDate)).setFullYear(s + a,
          0, 1), n[a] = angular.extend(this.createDateObject(t, this.formatYear), {
          uid: e.uniqueId + "-" + a
        });
        e.title = [n[0].label, n[i - 1].label].join(" - "), e.rows = this.split(n, r), e.columns = r
      }, this.compare = function(e, t) {
        return e.getFullYear() - t.getFullYear()
      }, this.handleKeyDown = function(e, t) {
        var n = this.activeDate.getFullYear();
        "left" === e ? n -= 1 : "up" === e ? n -= r : "right" === e ? n += 1 : "down" === e ? n += r : "pageup" === e || "pagedown" === e ? n += (
          "pageup" === e ? -1 : 1) * i : "home" === e ? n = o(this.activeDate.getFullYear()) : "end" === e && (n = o(this.activeDate
          .getFullYear()) + i - 1), this.activeDate.setFullYear(n)
      }
    }]).directive("uibDatepicker", (function() {
      return {
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/datepicker/datepicker.html"
        },
        scope: {
          datepickerOptions: "=?"
        },
        require: ["uibDatepicker", "^ngModel"],
        restrict: "A",
        controller: "UibDatepickerController",
        controllerAs: "datepicker",
        link: function(e, t, n, r) {
          var i = r[0],
            o = r[1];
          i.init(o)
        }
      }
    })).directive("uibDaypicker", (function() {
      return {
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/datepicker/day.html"
        },
        require: ["^uibDatepicker", "uibDaypicker"],
        restrict: "A",
        controller: "UibDaypickerController",
        link: function(e, t, n, r) {
          var i = r[0];
          r[1].init(i)
        }
      }
    })).directive("uibMonthpicker", (function() {
      return {
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/datepicker/month.html"
        },
        require: ["^uibDatepicker", "uibMonthpicker"],
        restrict: "A",
        controller: "UibMonthpickerController",
        link: function(e, t, n, r) {
          var i = r[0];
          r[1].init(i)
        }
      }
    })).directive("uibYearpicker", (function() {
      return {
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/datepicker/year.html"
        },
        require: ["^uibDatepicker", "uibYearpicker"],
        restrict: "A",
        controller: "UibYearpickerController",
        link: function(e, t, n, r) {
          var i = r[0];
          angular.extend(i, r[1]), i.yearpickerInit(), i.refreshView()
        }
      }
    })), angular.module("ui.bootstrap.position", []).factory("$uibPosition", ["$document", "$window", function(e, t) {
      var n, r, i = {
          normal: /(auto|scroll)/,
          hidden: /(auto|scroll|hidden)/
        },
        o = {
          auto: /\s?auto?\s?/i,
          primary: /^(top|bottom|left|right)$/,
          secondary: /^(top|bottom|left|right|center)$/,
          vertical: /^(top|bottom)$/
        },
        a = /(HTML|BODY)/;
      return {
        getRawNode: function(e) {
          return e.nodeName ? e : e[0] || e
        },
        parseStyle: function(e) {
          return e = parseFloat(e), isFinite(e) ? e : 0
        },
        offsetParent: function(n) {
          var r, i = (n = this.getRawNode(n)).offsetParent || e[0].documentElement;
          for (; i && i !== e[0].documentElement && (r = i, "static" === (t.getComputedStyle(r).position || "static"));) i = i.offsetParent;
          return i || e[0].documentElement
        },
        scrollbarWidth: function(i) {
          if (i) {
            if (angular.isUndefined(r)) {
              var o = e.find("body");
              o.addClass("uib-position-body-scrollbar-measure"), r = t.innerWidth - o[0].clientWidth, r = isFinite(r) ? r : 0, o.removeClass(
                "uib-position-body-scrollbar-measure")
            }
            return r
          }
          if (angular.isUndefined(n)) {
            var a = angular.element('<div class="uib-position-scrollbar-measure"></div>');
            e.find("body").append(a), n = a[0].offsetWidth - a[0].clientWidth, n = isFinite(n) ? n : 0, a.remove()
          }
          return n
        },
        scrollbarPadding: function(e) {
          e = this.getRawNode(e);
          var n = t.getComputedStyle(e),
            r = this.parseStyle(n.paddingRight),
            i = this.parseStyle(n.paddingBottom),
            o = this.scrollParent(e, !1, !0),
            s = this.scrollbarWidth(a.test(o.tagName));
          return {
            scrollbarWidth: s,
            widthOverflow: o.scrollWidth > o.clientWidth,
            right: r + s,
            originalRight: r,
            heightOverflow: o.scrollHeight > o.clientHeight,
            bottom: i + s,
            originalBottom: i
          }
        },
        isScrollable: function(e, n) {
          e = this.getRawNode(e);
          var r = n ? i.hidden : i.normal,
            o = t.getComputedStyle(e);
          return r.test(o.overflow + o.overflowY + o.overflowX)
        },
        scrollParent: function(n, r, o) {
          n = this.getRawNode(n);
          var a = r ? i.hidden : i.normal,
            s = e[0].documentElement,
            u = t.getComputedStyle(n);
          if (o && a.test(u.overflow + u.overflowY + u.overflowX)) return n;
          var l = "absolute" === u.position,
            c = n.parentElement || s;
          if (c === s || "fixed" === u.position) return s;
          for (; c.parentElement && c !== s;) {
            var d = t.getComputedStyle(c);
            if (l && "static" !== d.position && (l = !1), !l && a.test(d.overflow + d.overflowY + d.overflowX)) break;
            c = c.parentElement
          }
          return c
        },
        position: function(n, r) {
          n = this.getRawNode(n);
          var i = this.offset(n);
          if (r) {
            var o = t.getComputedStyle(n);
            i.top -= this.parseStyle(o.marginTop), i.left -= this.parseStyle(o.marginLeft)
          }
          var a = this.offsetParent(n),
            s = {
              top: 0,
              left: 0
            };
          return a !== e[0].documentElement && ((s = this.offset(a)).top += a.clientTop - a.scrollTop, s.left += a.clientLeft - a.scrollLeft), {
            width: Math.round(angular.isNumber(i.width) ? i.width : n.offsetWidth),
            height: Math.round(angular.isNumber(i.height) ? i.height : n.offsetHeight),
            top: Math.round(i.top - s.top),
            left: Math.round(i.left - s.left)
          }
        },
        offset: function(n) {
          var r = (n = this.getRawNode(n)).getBoundingClientRect();
          return {
            width: Math.round(angular.isNumber(r.width) ? r.width : n.offsetWidth),
            height: Math.round(angular.isNumber(r.height) ? r.height : n.offsetHeight),
            top: Math.round(r.top + (t.pageYOffset || e[0].documentElement.scrollTop)),
            left: Math.round(r.left + (t.pageXOffset || e[0].documentElement.scrollLeft))
          }
        },
        viewportOffset: function(n, r, i) {
          i = !1 !== i;
          var o = (n = this.getRawNode(n)).getBoundingClientRect(),
            a = {
              top: 0,
              left: 0,
              bottom: 0,
              right: 0
            },
            s = r ? e[0].documentElement : this.scrollParent(n),
            u = s.getBoundingClientRect();
          if (a.top = u.top + s.clientTop, a.left = u.left + s.clientLeft, s === e[0].documentElement && (a.top += t.pageYOffset, a.left += t
              .pageXOffset), a.bottom = a.top + s.clientHeight, a.right = a.left + s.clientWidth, i) {
            var l = t.getComputedStyle(s);
            a.top += this.parseStyle(l.paddingTop), a.bottom -= this.parseStyle(l.paddingBottom), a.left += this.parseStyle(l.paddingLeft), a
              .right -= this.parseStyle(l.paddingRight)
          }
          return {
            top: Math.round(o.top - a.top),
            bottom: Math.round(a.bottom - o.bottom),
            left: Math.round(o.left - a.left),
            right: Math.round(a.right - o.right)
          }
        },
        parsePlacement: function(e) {
          var t = o.auto.test(e);
          return t && (e = e.replace(o.auto, "")), (e = e.split("-"))[0] = e[0] || "top", o.primary.test(e[0]) || (e[0] = "top"), e[1] = e[1] ||
            "center", o.secondary.test(e[1]) || (e[1] = "center"), e[2] = !!t, e
        },
        positionElements: function(e, n, r, i) {
          e = this.getRawNode(e), n = this.getRawNode(n);
          var a = angular.isDefined(n.offsetWidth) ? n.offsetWidth : n.prop("offsetWidth"),
            s = angular.isDefined(n.offsetHeight) ? n.offsetHeight : n.prop("offsetHeight");
          r = this.parsePlacement(r);
          var u = i ? this.offset(e) : this.position(e),
            l = {
              top: 0,
              left: 0,
              placement: ""
            };
          if (r[2]) {
            var c = this.viewportOffset(e, i),
              d = t.getComputedStyle(n),
              p = {
                width: a + Math.round(Math.abs(this.parseStyle(d.marginLeft) + this.parseStyle(d.marginRight))),
                height: s + Math.round(Math.abs(this.parseStyle(d.marginTop) + this.parseStyle(d.marginBottom)))
              };
            if (r[0] = "top" === r[0] && p.height > c.top && p.height <= c.bottom ? "bottom" : "bottom" === r[0] && p.height > c.bottom && p
              .height <= c.top ? "top" : "left" === r[0] && p.width > c.left && p.width <= c.right ? "right" : "right" === r[0] && p.width > c
              .right && p.width <= c.left ? "left" : r[0], r[1] = "top" === r[1] && p.height - u.height > c.bottom && p.height - u.height <= c
              .top ? "bottom" : "bottom" === r[1] && p.height - u.height > c.top && p.height - u.height <= c.bottom ? "top" : "left" === r[1] && p
              .width - u.width > c.right && p.width - u.width <= c.left ? "right" : "right" === r[1] && p.width - u.width > c.left && p.width - u
              .width <= c.right ? "left" : r[1], "center" === r[1])
              if (o.vertical.test(r[0])) {
                var f = u.width / 2 - a / 2;
                c.left + f < 0 && p.width - u.width <= c.right ? r[1] = "left" : c.right + f < 0 && p.width - u.width <= c.left && (r[1] =
                  "right")
              } else {
                var h = u.height / 2 - p.height / 2;
                c.top + h < 0 && p.height - u.height <= c.bottom ? r[1] = "top" : c.bottom + h < 0 && p.height - u.height <= c.top && (r[1] =
                  "bottom")
              }
          }
          switch (r[0]) {
            case "top":
              l.top = u.top - s;
              break;
            case "bottom":
              l.top = u.top + u.height;
              break;
            case "left":
              l.left = u.left - a;
              break;
            case "right":
              l.left = u.left + u.width
          }
          switch (r[1]) {
            case "top":
              l.top = u.top;
              break;
            case "bottom":
              l.top = u.top + u.height - s;
              break;
            case "left":
              l.left = u.left;
              break;
            case "right":
              l.left = u.left + u.width - a;
              break;
            case "center":
              o.vertical.test(r[0]) ? l.left = u.left + u.width / 2 - a / 2 : l.top = u.top + u.height / 2 - s / 2
          }
          return l.top = Math.round(l.top), l.left = Math.round(l.left), l.placement = "center" === r[1] ? r[0] : r[0] + "-" + r[1], l
        },
        adjustTop: function(e, t, n, r) {
          if (-1 !== e.indexOf("top") && n !== r) return {
            top: t.top - r + "px"
          }
        },
        positionArrow: function(e, n) {
          var r = (e = this.getRawNode(e)).querySelector(".tooltip-inner, .popover-inner");
          if (r) {
            var i = angular.element(r).hasClass("tooltip-inner"),
              a = i ? e.querySelector(".tooltip-arrow") : e.querySelector(".arrow");
            if (a) {
              var s = {
                top: "",
                bottom: "",
                left: "",
                right: ""
              };
              if ("center" !== (n = this.parsePlacement(n))[1]) {
                var u = "border-" + n[0] + "-width",
                  l = t.getComputedStyle(a)[u],
                  c = "border-";
                o.vertical.test(n[0]) ? c += n[0] + "-" + n[1] : c += n[1] + "-" + n[0], c += "-radius";
                var d = t.getComputedStyle(i ? r : e)[c];
                switch (n[0]) {
                  case "top":
                    s.bottom = i ? "0" : "-" + l;
                    break;
                  case "bottom":
                    s.top = i ? "0" : "-" + l;
                    break;
                  case "left":
                    s.right = i ? "0" : "-" + l;
                    break;
                  case "right":
                    s.left = i ? "0" : "-" + l
                }
                s[n[1]] = d, angular.element(a).css(s)
              } else angular.element(a).css(s)
            }
          }
        }
      }
    }]), angular.module("ui.bootstrap.datepickerPopup", ["ui.bootstrap.datepicker", "ui.bootstrap.position"]).value("$datepickerPopupLiteralWarning", !
      0).constant("uibDatepickerPopupConfig", {
      altInputFormats: [],
      appendToBody: !1,
      clearText: "Clear",
      closeOnDateSelection: !0,
      closeText: "Done",
      currentText: "Today",
      datepickerPopup: "yyyy-MM-dd",
      datepickerPopupTemplateUrl: "uib/template/datepickerPopup/popup.html",
      datepickerTemplateUrl: "uib/template/datepicker/datepicker.html",
      html5Types: {
        date: "yyyy-MM-dd",
        "datetime-local": "yyyy-MM-ddTHH:mm:ss.sss",
        month: "yyyy-MM"
      },
      onOpenFocus: !0,
      showButtonBar: !0,
      placement: "auto bottom-left"
    }).controller("UibDatepickerPopupController", ["$scope", "$element", "$attrs", "$compile", "$log", "$parse", "$window", "$document", "$rootScope",
      "$uibPosition", "dateFilter", "uibDateParser", "uibDatepickerPopupConfig", "$timeout", "uibDatepickerConfig", "$datepickerPopupLiteralWarning",
      function(e, t, n, r, i, o, a, s, u, l, c, d, p, f, h, g) {
        var m, v, $, b, y, w, x, C, S, A, T, k, D, M = !1,
          E = [];

        function O(t) {
          var n = d.parse(t, m, e.date);
          if (isNaN(n))
            for (var r = 0; r < D.length; r++)
              if (n = d.parse(t, D[r], e.date), !isNaN(n)) return n;
          return n
        }

        function N(e) {
          if (angular.isNumber(e) && (e = new Date(e)), !e) return null;
          if (angular.isDate(e) && !isNaN(e)) return e;
          if (angular.isString(e)) {
            var t = O(e);
            if (!isNaN(t)) return d.toTimezone(t, T.getOption("timezone"))
          }
          return T.getOption("allowInvalid") ? e : void 0
        }

        function P(e, t) {
          var r = e || t;
          return !n.ngRequired && !r || (angular.isNumber(r) && (r = new Date(r)), !r || (!(!angular.isDate(r) || isNaN(r)) || !!angular.isString(
            r) && !isNaN(O(r))))
        }

        function q(n) {
          if (e.isOpen || !e.disabled) {
            var r = k[0],
              i = t[0].contains(n.target),
              o = void 0 !== r.contains && r.contains(n.target);
            !e.isOpen || i || o || e.$apply((function() {
              e.isOpen = !1
            }))
          }
        }

        function I(n) {
          27 === n.which && e.isOpen ? (n.preventDefault(), n.stopPropagation(), e.$apply((function() {
            e.isOpen = !1
          })), t[0].focus()) : 40 !== n.which || e.isOpen || (n.preventDefault(), n.stopPropagation(), e.$apply((function() {
            e.isOpen = !0
          })))
        }

        function L() {
          if (e.isOpen) {
            var r = angular.element(k[0].querySelector(".uib-datepicker-popup")),
              i = n.popupPlacement ? n.popupPlacement : p.placement,
              o = l.positionElements(t, r, i, $);
            r.css({
              top: o.top + "px",
              left: o.left + "px"
            }), r.hasClass("uib-position-measure") && r.removeClass("uib-position-measure")
          }
        }
        this.init = function(i) {
          if (T = function(e) {
              var t;
              angular.version.minor < 6 ? (t = angular.isObject(e.$options) ? e.$options : {
                timezone: null
              }).getOption = function(e) {
                return t[e]
              } : t = e.$options;
              return t
            }(A = i), v = angular.isDefined(n.closeOnDateSelection) ? e.$parent.$eval(n.closeOnDateSelection) : p.closeOnDateSelection, $ = angular
            .isDefined(n.datepickerAppendToBody) ? e.$parent.$eval(n.datepickerAppendToBody) : p.appendToBody, b = angular.isDefined(n
            .onOpenFocus) ? e.$parent.$eval(n.onOpenFocus) : p.onOpenFocus, y = angular.isDefined(n.datepickerPopupTemplateUrl) ? n
            .datepickerPopupTemplateUrl : p.datepickerPopupTemplateUrl, w = angular.isDefined(n.datepickerTemplateUrl) ? n.datepickerTemplateUrl : p
            .datepickerTemplateUrl, D = angular.isDefined(n.altInputFormats) ? e.$parent.$eval(n.altInputFormats) : p.altInputFormats, e
            .showButtonBar = angular.isDefined(n.showButtonBar) ? e.$parent.$eval(n.showButtonBar) : p.showButtonBar, p.html5Types[n.type] ? (m = p
              .html5Types[n.type], M = !0) : (m = n.uibDatepickerPopup || p.datepickerPopup, n.$observe("uibDatepickerPopup", (function(e, t) {
              var n = e || p.datepickerPopup;
              if (n !== m && (m = n, A.$modelValue = null, !m)) throw new Error("uibDatepickerPopup must have a date format specified.")
            }))), !m) throw new Error("uibDatepickerPopup must have a date format specified.");
          if (M && n.uibDatepickerPopup) throw new Error("HTML5 date input types do not support custom formats.");
          (x = angular.element("<div uib-datepicker-popup-wrap><div uib-datepicker></div></div>")).attr({
              "ng-model": "date",
              "ng-change": "dateSelection(date)",
              "template-url": y
            }), (C = angular.element(x.children()[0])).attr("template-url", w), e.datepickerOptions || (e.datepickerOptions = {}), M && "month" ===
            n.type && (e.datepickerOptions.datepickerMode = "month", e.datepickerOptions.minMode = "month"), C.attr("datepicker-options",
              "datepickerOptions"), M ? A.$formatters.push((function(t) {
              return e.date = d.fromTimezone(t, T.getOption("timezone")), t
            })) : (A.$$parserName = "date", A.$validators.date = P, A.$parsers.unshift(N), A.$formatters.push((function(t) {
              return A.$isEmpty(t) ? (e.date = t, t) : (angular.isNumber(t) && (t = new Date(t)), e.date = d.fromTimezone(t, T.getOption(
                "timezone")), d.filter(e.date, m))
            }))), A.$viewChangeListeners.push((function() {
              e.date = O(A.$viewValue)
            })), t.on("keydown", I), k = r(x)(e), x.remove(), $ ? s.find("body").append(k) : t.after(k), e.$on("$destroy", (function() {
              for (!0 === e.isOpen && (u.$$phase || e.$apply((function() {
                  e.isOpen = !1
                }))), k.remove(), t.off("keydown", I), s.off("click", q), S && S.off("scroll", L), angular.element(a).off("resize", L); E
                .length;) E.shift()()
            }))
        }, e.getText = function(t) {
          return e[t + "Text"] || p[t + "Text"]
        }, e.isDisabled = function(t) {
          "today" === t && (t = d.fromTimezone(new Date, T.getOption("timezone")));
          var n = {};
          return angular.forEach(["minDate", "maxDate"], (function(t) {
            e.datepickerOptions[t] ? angular.isDate(e.datepickerOptions[t]) ? n[t] = new Date(e.datepickerOptions[t]) : (g && i.warn(
              "Literal date support has been deprecated, please switch to date object usage"), n[t] = new Date(c(e.datepickerOptions[t],
              "medium"))) : n[t] = null
          })), e.datepickerOptions && n.minDate && e.compare(t, n.minDate) < 0 || n.maxDate && e.compare(t, n.maxDate) > 0
        }, e.compare = function(e, t) {
          return new Date(e.getFullYear(), e.getMonth(), e.getDate()) - new Date(t.getFullYear(), t.getMonth(), t.getDate())
        }, e.dateSelection = function(n) {
          e.date = n;
          var r = e.date ? d.filter(e.date, m) : null;
          t.val(r), A.$setViewValue(r), v && (e.isOpen = !1, t[0].focus())
        }, e.keydown = function(n) {
          27 === n.which && (n.stopPropagation(), e.isOpen = !1, t[0].focus())
        }, e.select = function(t, n) {
          if (n.stopPropagation(), "today" === t) {
            var r = new Date;
            angular.isDate(e.date) ? (t = new Date(e.date)).setFullYear(r.getFullYear(), r.getMonth(), r.getDate()) : (t = d.fromTimezone(r, T
              .getOption("timezone"))).setHours(0, 0, 0, 0)
          }
          e.dateSelection(t)
        }, e.close = function(n) {
          n.stopPropagation(), e.isOpen = !1, t[0].focus()
        }, e.disabled = angular.isDefined(n.disabled) || !1, n.ngDisabled && E.push(e.$parent.$watch(o(n.ngDisabled), (function(t) {
          e.disabled = t
        }))), e.$watch("isOpen", (function(r) {
          r ? e.disabled ? e.isOpen = !1 : f((function() {
            L(), b && e.$broadcast("uib:datepicker.focus"), s.on("click", q);
            var r = n.popupPlacement ? n.popupPlacement : p.placement;
            $ || l.parsePlacement(r)[2] ? (S = S || angular.element(l.scrollParent(t))) && S.on("scroll", L) : S = null, angular.element(
              a).on("resize", L)
          }), 0, !1) : (s.off("click", q), S && S.off("scroll", L), angular.element(a).off("resize", L))
        })), e.$on("uib:datepicker.mode", (function() {
          f(L, 0, !1)
        }))
      }
    ]).directive("uibDatepickerPopup", (function() {
      return {
        require: ["ngModel", "uibDatepickerPopup"],
        controller: "UibDatepickerPopupController",
        scope: {
          datepickerOptions: "=?",
          isOpen: "=?",
          currentText: "@",
          clearText: "@",
          closeText: "@"
        },
        link: function(e, t, n, r) {
          var i = r[0];
          r[1].init(i)
        }
      }
    })).directive("uibDatepickerPopupWrap", (function() {
      return {
        restrict: "A",
        transclude: !0,
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/datepickerPopup/popup.html"
        }
      }
    })), angular.module("ui.bootstrap.debounce", []).factory("$$debounce", ["$timeout", function(e) {
      return function(t, n) {
        var r;
        return function() {
          var i = this,
            o = Array.prototype.slice.call(arguments);
          r && e.cancel(r), r = e((function() {
            t.apply(i, o)
          }), n)
        }
      }
    }]), angular.module("ui.bootstrap.multiMap", []).factory("$$multiMap", (function() {
      return {
        createNew: function() {
          var e = {};
          return {
            entries: function() {
              return Object.keys(e).map((function(t) {
                return {
                  key: t,
                  value: e[t]
                }
              }))
            },
            get: function(t) {
              return e[t]
            },
            hasKey: function(t) {
              return !!e[t]
            },
            keys: function() {
              return Object.keys(e)
            },
            put: function(t, n) {
              e[t] || (e[t] = []), e[t].push(n)
            },
            remove: function(t, n) {
              var r = e[t];
              if (r) {
                var i = r.indexOf(n); - 1 !== i && r.splice(i, 1), r.length || delete e[t]
              }
            }
          }
        }
      }
    })), angular.module("ui.bootstrap.dropdown", ["ui.bootstrap.multiMap", "ui.bootstrap.position"]).constant("uibDropdownConfig", {
      appendToOpenClass: "uib-dropdown-open",
      openClass: "open"
    }).service("uibDropdownService", ["$document", "$rootScope", "$$multiMap", function(e, t, n) {
      var r = null,
        i = n.createNew();
      this.isOnlyOpen = function(e, t) {
        var n = i.get(t);
        if (n && n.reduce((function(t, n) {
            return n.scope === e ? n : t
          }), {})) return 1 === n.length;
        return !1
      }, this.open = function(t, n, a) {
        if (r || e.on("click", o), r && r !== t && (r.isOpen = !1), r = t, a) {
          var s = i.get(a);
          if (s) - 1 === s.map((function(e) {
            return e.scope
          })).indexOf(t) && i.put(a, {
            scope: t
          });
          else i.put(a, {
            scope: t
          })
        }
      }, this.close = function(t, n, a) {
        if (r === t && (e.off("click", o), e.off("keydown", this.keybindFilter), r = null), a) {
          var s = i.get(a);
          if (s) {
            var u = s.reduce((function(e, n) {
              return n.scope === t ? n : e
            }), {});
            u && i.remove(a, u)
          }
        }
      };
      var o = function(e) {
        if (r && r.isOpen && !(e && "disabled" === r.getAutoClose() || e && 3 === e.which)) {
          var n = r.getToggleElement();
          if (!(e && n && n[0].contains(e.target))) {
            var i = r.getDropdownElement();
            e && "outsideClick" === r.getAutoClose() && i && i[0].contains(e.target) || (r.focusToggleElement(), r.isOpen = !1, t.$$phase || r
              .$apply())
          }
        }
      };
      this.keybindFilter = function(e) {
        if (r) {
          var t = r.getDropdownElement(),
            n = r.getToggleElement(),
            i = t && t[0].contains(e.target),
            a = n && n[0].contains(e.target);
          27 === e.which ? (e.stopPropagation(), r.focusToggleElement(), o()) : r.isKeynavEnabled() && -1 !== [38, 40].indexOf(e.which) && r
            .isOpen && (i || a) && (e.preventDefault(), e.stopPropagation(), r.focusDropdownEntry(e.which))
        }
      }
    }]).controller("UibDropdownController", ["$scope", "$element", "$attrs", "$parse", "uibDropdownConfig", "uibDropdownService", "$animate",
      "$uibPosition", "$document", "$compile", "$templateRequest",
      function(e, t, n, r, i, o, a, s, u, l, c) {
        var d, p, f = this,
          h = e.$new(),
          g = i.appendToOpenClass,
          m = i.openClass,
          v = angular.noop,
          $ = n.onToggle ? r(n.onToggle) : angular.noop,
          b = !1,
          y = u.find("body");

        function w() {
          t.append(f.dropdownMenu)
        }
        t.addClass("dropdown"), this.init = function() {
          n.isOpen && (p = r(n.isOpen), v = p.assign, e.$watch(p, (function(e) {
            h.isOpen = !!e
          }))), b = angular.isDefined(n.keyboardNav)
        }, this.toggle = function(e) {
          return h.isOpen = arguments.length ? !!e : !h.isOpen, angular.isFunction(v) && v(h, h.isOpen), h.isOpen
        }, this.isOpen = function() {
          return h.isOpen
        }, h.getToggleElement = function() {
          return f.toggleElement
        }, h.getAutoClose = function() {
          return n.autoClose || "always"
        }, h.getElement = function() {
          return t
        }, h.isKeynavEnabled = function() {
          return b
        }, h.focusDropdownEntry = function(e) {
          var n = f.dropdownMenu ? angular.element(f.dropdownMenu).find("a") : t.find("ul").eq(0).find("a");
          switch (e) {
            case 40:
              angular.isNumber(f.selectedOption) ? f.selectedOption = f.selectedOption === n.length - 1 ? f.selectedOption : f.selectedOption + 1 :
                f.selectedOption = 0;
              break;
            case 38:
              angular.isNumber(f.selectedOption) ? f.selectedOption = 0 === f.selectedOption ? 0 : f.selectedOption - 1 : f.selectedOption = n
                .length - 1
          }
          n[f.selectedOption].focus()
        }, h.getDropdownElement = function() {
          return f.dropdownMenu
        }, h.focusToggleElement = function() {
          f.toggleElement && f.toggleElement[0].focus()
        }, h.$watch("isOpen", (function(i, p) {
          var b = null,
            x = !1;
          if (angular.isDefined(n.dropdownAppendTo)) {
            var C = r(n.dropdownAppendTo)(h);
            C && (b = angular.element(C))
          }
          angular.isDefined(n.dropdownAppendToBody) && (!1 !== r(n.dropdownAppendToBody)(h) && (x = !0));
          if (x && !b && (b = y), b && f.dropdownMenu && (i ? (b.append(f.dropdownMenu), t.on("$destroy", w)) : (t.off("$destroy", w), w())),
            b && f.dropdownMenu) {
            var S, A, T, k = s.positionElements(t, f.dropdownMenu, "bottom-left", !0),
              D = 0;
            if (S = {
                top: k.top + "px",
                display: i ? "block" : "none"
              }, (A = f.dropdownMenu.hasClass("dropdown-menu-right")) ? (S.left = "auto", (T = s.scrollbarPadding(b)).heightOverflow && T
                .scrollbarWidth && (D = T.scrollbarWidth), S.right = window.innerWidth - D - (k.left + t.prop("offsetWidth")) + "px") : (S
                .left = k.left + "px", S.right = "auto"), !x) {
              var M = s.offset(b);
              S.top = k.top - M.top + "px", A ? S.right = window.innerWidth - (k.left - M.left + t.prop("offsetWidth")) + "px" : S.left = k
                .left - M.left + "px"
            }
            f.dropdownMenu.css(S)
          }
          var E = b || t,
            O = b ? g : m,
            N = E.hasClass(O),
            P = o.isOnlyOpen(e, b);
          N === !i && a[b ? P ? "removeClass" : "addClass" : i ? "addClass" : "removeClass"](E, O).then((function() {
            angular.isDefined(i) && i !== p && $(e, {
              open: !!i
            })
          }));
          if (i) f.dropdownMenuTemplateUrl ? c(f.dropdownMenuTemplateUrl).then((function(e) {
            d = h.$new(), l(e.trim())(d, (function(e) {
              var t = e;
              f.dropdownMenu.replaceWith(t), f.dropdownMenu = t, u.on("keydown", o.keybindFilter)
            }))
          })) : u.on("keydown", o.keybindFilter), h.focusToggleElement(), o.open(h, t, b);
          else {
            if (o.close(h, t, b), f.dropdownMenuTemplateUrl) {
              d && d.$destroy();
              var q = angular.element('<ul class="dropdown-menu"></ul>');
              f.dropdownMenu.replaceWith(q), f.dropdownMenu = q
            }
            f.selectedOption = null
          }
          angular.isFunction(v) && v(e, i)
        }))
      }
    ]).directive("uibDropdown", (function() {
      return {
        controller: "UibDropdownController",
        link: function(e, t, n, r) {
          r.init()
        }
      }
    })).directive("uibDropdownMenu", (function() {
      return {
        restrict: "A",
        require: "?^uibDropdown",
        link: function(e, t, n, r) {
          if (r && !angular.isDefined(n.dropdownNested)) {
            t.addClass("dropdown-menu");
            var i = n.templateUrl;
            i && (r.dropdownMenuTemplateUrl = i), r.dropdownMenu || (r.dropdownMenu = t)
          }
        }
      }
    })).directive("uibDropdownToggle", (function() {
      return {
        require: "?^uibDropdown",
        link: function(e, t, n, r) {
          if (r) {
            t.addClass("dropdown-toggle"), r.toggleElement = t;
            var i = function(i) {
              i.preventDefault(), t.hasClass("disabled") || n.disabled || e.$apply((function() {
                r.toggle()
              }))
            };
            t.on("click", i), t.attr({
              "aria-haspopup": !0,
              "aria-expanded": !1
            }), e.$watch(r.isOpen, (function(e) {
              t.attr("aria-expanded", !!e)
            })), e.$on("$destroy", (function() {
              t.off("click", i)
            }))
          }
        }
      }
    })), angular.module("ui.bootstrap.stackedMap", []).factory("$$stackedMap", (function() {
      return {
        createNew: function() {
          var e = [];
          return {
            add: function(t, n) {
              e.push({
                key: t,
                value: n
              })
            },
            get: function(t) {
              for (var n = 0; n < e.length; n++)
                if (t === e[n].key) return e[n]
            },
            keys: function() {
              for (var t = [], n = 0; n < e.length; n++) t.push(e[n].key);
              return t
            },
            top: function() {
              return e[e.length - 1]
            },
            remove: function(t) {
              for (var n = -1, r = 0; r < e.length; r++)
                if (t === e[r].key) {
                  n = r;
                  break
                } return e.splice(n, 1)[0]
            },
            removeTop: function() {
              return e.pop()
            },
            length: function() {
              return e.length
            }
          }
        }
      }
    })), angular.module("ui.bootstrap.modal", ["ui.bootstrap.multiMap", "ui.bootstrap.stackedMap", "ui.bootstrap.position"]).provider("$uibResolve", (
      function() {
        var e = this;
        this.resolver = null, this.setResolver = function(e) {
          this.resolver = e
        }, this.$get = ["$injector", "$q", function(t, n) {
          var r = e.resolver ? t.get(e.resolver) : null;
          return {
            resolve: function(e, i, o, a) {
              if (r) return r.resolve(e, i, o, a);
              var s = [];
              return angular.forEach(e, (function(e) {
                angular.isFunction(e) || angular.isArray(e) ? s.push(n.resolve(t.invoke(e))) : angular.isString(e) ? s.push(n.resolve(t
                  .get(e))) : s.push(n.resolve(e))
              })), n.all(s).then((function(t) {
                var n = {},
                  r = 0;
                return angular.forEach(e, (function(e, i) {
                  n[i] = t[r++]
                })), n
              }))
            }
          }
        }]
      })).directive("uibModalBackdrop", ["$animate", "$injector", "$uibModalStack", function(e, t, n) {
      return {
        restrict: "A",
        compile: function(e, t) {
          return e.addClass(t.backdropClass), r
        }
      };

      function r(t, r, i) {
        i.modalInClass && (e.addClass(r, i.modalInClass), t.$on(n.NOW_CLOSING_EVENT, (function(n, o) {
          var a = o();
          t.modalOptions.animation ? e.removeClass(r, i.modalInClass).then(a) : a()
        })))
      }
    }]).directive("uibModalWindow", ["$uibModalStack", "$q", "$animateCss", "$document", function(e, t, n, r) {
      return {
        scope: {
          index: "@"
        },
        restrict: "A",
        transclude: !0,
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/modal/window.html"
        },
        link: function(i, o, a) {
          o.addClass(a.windowTopClass || ""), i.size = a.size, i.close = function(t) {
            var n = e.getTop();
            n && n.value.backdrop && "static" !== n.value.backdrop && t.target === t.currentTarget && (t.preventDefault(), t.stopPropagation(),
              e.dismiss(n.key, "backdrop click"))
          }, o.on("click", i.close), i.$isRendered = !0;
          var s = t.defer();
          i.$$postDigest((function() {
            s.resolve()
          })), s.promise.then((function() {
            var s = null;
            a.modalInClass && (s = n(o, {
              addClass: a.modalInClass
            }).start(), i.$on(e.NOW_CLOSING_EVENT, (function(e, t) {
              var r = t();
              n(o, {
                removeClass: a.modalInClass
              }).start().then(r)
            }))), t.when(s).then((function() {
              var t = e.getTop();
              if (t && e.modalRendered(t.key), !r[0].activeElement || !o[0].contains(r[0].activeElement)) {
                var n = o[0].querySelector("[autofocus]");
                n ? n.focus() : o[0].focus()
              }
            }))
          }))
        }
      }
    }]).directive("uibModalAnimationClass", (function() {
      return {
        compile: function(e, t) {
          t.modalAnimation && e.addClass(t.uibModalAnimationClass)
        }
      }
    })).directive("uibModalTransclude", ["$animate", function(e) {
      return {
        link: function(t, n, r, i, o) {
          o(t.$parent, (function(t) {
            n.empty(), e.enter(t, n)
          }))
        }
      }
    }]).factory("$uibModalStack", ["$animate", "$animateCss", "$document", "$compile", "$rootScope", "$q", "$$multiMap", "$$stackedMap", "$uibPosition",
      function(e, t, n, r, i, o, a, s, u) {
        var l, c, d, p = "modal-open",
          f = s.createNew(),
          h = a.createNew(),
          g = {
            NOW_CLOSING_EVENT: "modal.stack.now-closing"
          },
          m = 0,
          v = null,
          $ = "data-bootstrap-modal-aria-hidden-count",
          b = /[A-Z]/g;

        function y() {
          for (var e = -1, t = f.keys(), n = 0; n < t.length; n++) f.get(t[n]).value.backdrop && (e = n);
          return e > -1 && e < m && (e = m), e
        }

        function w(e, t) {
          var n = f.get(e).value,
            r = n.appendTo;
          f.remove(e), (v = f.top()) && (m = parseInt(v.value.modalDomEl.attr("index"), 10)), C(n.modalDomEl, n.modalScope, (function() {
              var t = n.openedClass || p;
              h.remove(t, e);
              var i = h.hasKey(t);
              r.toggleClass(t, i), !i && d && d.heightOverflow && d.scrollbarWidth && (d.originalRight ? r.css({
                paddingRight: d.originalRight + "px"
              }) : r.css({
                paddingRight: ""
              }), d = null), x(!0)
            }), n.closedDeferred),
            function() {
              if (l && -1 === y()) {
                C(l, c, (function() {
                  null
                })), l = void 0, c = void 0
              }
            }(), t && t.focus ? t.focus() : r.focus && r.focus()
        }

        function x(e) {
          var t;
          f.length() > 0 && (t = f.top().value).modalDomEl.toggleClass(t.windowTopClass || "", e)
        }

        function C(t, n, r, i) {
          var a, s = null;
          return n.$broadcast(g.NOW_CLOSING_EVENT, (function() {
            return a || (a = o.defer(), s = a.promise),
              function() {
                a.resolve()
              }
          })), o.when(s).then((function o() {
            if (o.done) return;
            o.done = !0, e.leave(t).then((function() {
              r && r(), t.remove(), i && i.resolve()
            })), n.$destroy()
          }))
        }

        function S(e) {
          if (e.isDefaultPrevented()) return e;
          var t = f.top();
          if (t) switch (e.which) {
            case 27:
              t.value.keyboard && (e.preventDefault(), i.$apply((function() {
                g.dismiss(t.key, "escape key press")
              })));
              break;
            case 9:
              var n = g.loadFocusElementList(t),
                r = !1;
              e.shiftKey ? (g.isFocusInFirstItem(e, n) || g.isModalFocused(e, t)) && (r = g.focusLastFocusableElement(n)) : g.isFocusInLastItem(e,
                n) && (r = g.focusFirstFocusableElement(n)), r && (e.preventDefault(), e.stopPropagation())
          }
        }

        function A(e, t, n) {
          return !e.value.modalScope.$broadcast("modal.closing", t, n).defaultPrevented
        }

        function T() {
          Array.prototype.forEach.call(document.querySelectorAll("[" + $ + "]"), (function(e) {
            var t = parseInt(e.getAttribute($), 10) - 1;
            e.setAttribute($, t), t || (e.removeAttribute($), e.removeAttribute("aria-hidden"))
          }))
        }
        return i.$watch(y, (function(e) {
          c && (c.index = e)
        })), n.on("keydown", S), i.$on("$destroy", (function() {
          n.off("keydown", S)
        })), g.open = function(t, o) {
          var a = n[0].activeElement,
            s = o.openedClass || p;
          x(!1), v = f.top(), f.add(t, {
            deferred: o.deferred,
            renderDeferred: o.renderDeferred,
            closedDeferred: o.closedDeferred,
            modalScope: o.scope,
            backdrop: o.backdrop,
            keyboard: o.keyboard,
            openedClass: o.openedClass,
            windowTopClass: o.windowTopClass,
            animation: o.animation,
            appendTo: o.appendTo
          }), h.put(s, t);
          var g, w = o.appendTo,
            C = y();
          C >= 0 && !l && ((c = i.$new(!0)).modalOptions = o, c.index = C, (l = angular.element('<div uib-modal-backdrop="modal-backdrop"></div>'))
            .attr({
              class: "modal-backdrop",
              "ng-style": "{'z-index': 1040 + (index && 1 || 0) + index*10}",
              "uib-modal-animation-class": "fade",
              "modal-in-class": "in"
            }), o.backdropClass && l.addClass(o.backdropClass), o.animation && l.attr("modal-animation", "true"), r(l)(c), e.enter(l, w), u
            .isScrollable(w) && (d = u.scrollbarPadding(w)).heightOverflow && d.scrollbarWidth && w.css({
              paddingRight: d.right + "px"
            })), o.component ? (g = document.createElement(o.component.name.replace(b, (function(e, t) {
            return (t ? "-" : "") + e.toLowerCase()
          }))), (g = angular.element(g)).attr({
            resolve: "$resolve",
            "modal-instance": "$uibModalInstance",
            close: "$close($value)",
            dismiss: "$dismiss($value)"
          })) : g = o.content, m = v ? parseInt(v.value.modalDomEl.attr("index"), 10) + 1 : 0;
          var S = angular.element('<div uib-modal-window="modal-window"></div>');
          S.attr({
              class: "modal",
              "template-url": o.windowTemplateUrl,
              "window-top-class": o.windowTopClass,
              role: "dialog",
              "aria-labelledby": o.ariaLabelledBy,
              "aria-describedby": o.ariaDescribedBy,
              size: o.size,
              index: m,
              animate: "animate",
              "ng-style": "{'z-index': 1050 + $$topModalIndex*10, display: 'block'}",
              tabindex: -1,
              "uib-modal-animation-class": "fade",
              "modal-in-class": "in"
            }).append(g), o.windowClass && S.addClass(o.windowClass), o.animation && S.attr("modal-animation", "true"), w.addClass(s), o.scope && (o
              .scope.$$topModalIndex = m), e.enter(r(S)(o.scope), w), f.top().value.modalDomEl = S, f.top().value.modalOpener = a,
            function e(t) {
              if (!t || "BODY" === t[0].tagName) return;
              return function(e) {
                var t = e.parent() ? e.parent().children() : [];
                return Array.prototype.filter.call(t, (function(t) {
                  return t !== e[0]
                }))
              }(t).forEach((function(e) {
                var t = "true" === e.getAttribute("aria-hidden"),
                  n = parseInt(e.getAttribute($), 10);
                n || (n = t ? 1 : 0), e.setAttribute($, n + 1), e.setAttribute("aria-hidden", "true")
              })), e(t.parent())
            }(S)
        }, g.close = function(e, t) {
          var n = f.get(e);
          return T(), n && A(n, t, !0) ? (n.value.modalScope.$$uibDestructionScheduled = !0, n.value.deferred.resolve(t), w(e, n.value.modalOpener),
            !0) : !n
        }, g.dismiss = function(e, t) {
          var n = f.get(e);
          return T(), n && A(n, t, !1) ? (n.value.modalScope.$$uibDestructionScheduled = !0, n.value.deferred.reject(t), w(e, n.value.modalOpener),
            !0) : !n
        }, g.dismissAll = function(e) {
          for (var t = this.getTop(); t && this.dismiss(t.key, e);) t = this.getTop()
        }, g.getTop = function() {
          return f.top()
        }, g.modalRendered = function(e) {
          var t = f.get(e);
          t && t.value.renderDeferred.resolve()
        }, g.focusFirstFocusableElement = function(e) {
          return e.length > 0 && (e[0].focus(), !0)
        }, g.focusLastFocusableElement = function(e) {
          return e.length > 0 && (e[e.length - 1].focus(), !0)
        }, g.isModalFocused = function(e, t) {
          if (e && t) {
            var n = t.value.modalDomEl;
            if (n && n.length) return (e.target || e.srcElement) === n[0]
          }
          return !1
        }, g.isFocusInFirstItem = function(e, t) {
          return t.length > 0 && (e.target || e.srcElement) === t[0]
        }, g.isFocusInLastItem = function(e, t) {
          return t.length > 0 && (e.target || e.srcElement) === t[t.length - 1]
        }, g.loadFocusElementList = function(e) {
          if (e) {
            var t = e.value.modalDomEl;
            if (t && t.length) {
              var n = t[0].querySelectorAll(
                "a[href], area[href], input:not([disabled]):not([tabindex='-1']), button:not([disabled]):not([tabindex='-1']),select:not([disabled]):not([tabindex='-1']), textarea:not([disabled]):not([tabindex='-1']), iframe, object, embed, *[tabindex]:not([tabindex='-1']), *[contenteditable=true]"
                );
              return n ? Array.prototype.filter.call(n, (function(e) {
                return function(e) {
                  return !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length)
                }(e)
              })) : n
            }
          }
        }, g
      }
    ]).provider("$uibModal", (function() {
      var e = {
        options: {
          animation: !0,
          backdrop: !0,
          keyboard: !0
        },
        $get: ["$rootScope", "$q", "$document", "$templateRequest", "$controller", "$uibResolve", "$uibModalStack", function(t, n, r, i, o, a,
        s) {
          var u = {};
          var l = null;
          return u.getPromiseChain = function() {
            return l
          }, u.open = function(u) {
            var c, d, p, f = n.defer(),
              h = n.defer(),
              g = n.defer(),
              m = n.defer(),
              v = {
                result: f.promise,
                opened: h.promise,
                closed: g.promise,
                rendered: m.promise,
                close: function(e) {
                  return s.close(v, e)
                },
                dismiss: function(e) {
                  return s.dismiss(v, e)
                }
              };
            if ((u = angular.extend({}, e.options, u)).resolve = u.resolve || {}, u.appendTo = u.appendTo || r.find("body").eq(0), !u
              .appendTo.length) throw new Error("appendTo element not found. Make sure that the element passed is in DOM.");
            if (!u.component && !u.template && !u.templateUrl) throw new Error(
              "One of component or template or templateUrl options is required.");

            function $() {
              return c
            }
            return c = u.component ? n.when(a.resolve(u.resolve, {}, null, null)) : n.all([(d = u, d.template ? n.when(d.template) : i(
                angular.isFunction(d.templateUrl) ? d.templateUrl() : d.templateUrl)), a.resolve(u.resolve, {}, null, null)]), p = l = n
              .all([l]).then($, $).then((function(e) {
                var n = u.scope || t,
                  r = n.$new();
                r.$close = v.close, r.$dismiss = v.dismiss, r.$on("$destroy", (function() {
                  r.$$uibDestructionScheduled || r.$dismiss("$uibUnscheduledDestruction")
                }));
                var i, a, l = {
                    scope: r,
                    deferred: f,
                    renderDeferred: m,
                    closedDeferred: g,
                    animation: u.animation,
                    backdrop: u.backdrop,
                    keyboard: u.keyboard,
                    backdropClass: u.backdropClass,
                    windowTopClass: u.windowTopClass,
                    windowClass: u.windowClass,
                    windowTemplateUrl: u.windowTemplateUrl,
                    ariaLabelledBy: u.ariaLabelledBy,
                    ariaDescribedBy: u.ariaDescribedBy,
                    size: u.size,
                    openedClass: u.openedClass,
                    appendTo: u.appendTo
                  },
                  c = {},
                  d = {};

                function p(t, n, i, o) {
                  t.$scope = r, t.$scope.$resolve = {}, i ? t.$scope.$uibModalInstance = v : t.$uibModalInstance = v;
                  var a = n ? e[1] : e;
                  angular.forEach(a, (function(e, n) {
                    o && (t[n] = e), t.$scope.$resolve[n] = e
                  }))
                }
                u.component ? (p(c, !1, !0, !1), c.name = u.component, l.component = c) : u.controller && (p(d, !0, !1, !0), a = o(u
                    .controller, d, !0, u.controllerAs), u.controllerAs && u.bindToController && ((i = a.instance).$close = r.$close,
                    i.$dismiss = r.$dismiss, angular.extend(i, {
                      $resolve: d.$scope.$resolve
                    }, n)), i = a(), angular.isFunction(i.$onInit) && i.$onInit()), u.component || (l.content = e[0]), s.open(v, l), h
                  .resolve(!0)
              }), (function(e) {
                h.reject(e), f.reject(e)
              })).finally((function() {
                l === p && (l = null)
              })), v
          }, u
        }]
      };
      return e
    })), angular.module("ui.bootstrap.paging", []).factory("uibPaging", ["$parse", function(e) {
      return {
        create: function(t, n, r) {
          t.setNumPages = r.numPages ? e(r.numPages).assign : angular.noop, t.ngModelCtrl = {
            $setViewValue: angular.noop
          }, t._watchers = [], t.init = function(e, i) {
            t.ngModelCtrl = e, t.config = i, e.$render = function() {
              t.render()
            }, r.itemsPerPage ? t._watchers.push(n.$parent.$watch(r.itemsPerPage, (function(e) {
              t.itemsPerPage = parseInt(e, 10), n.totalPages = t.calculateTotalPages(), t.updatePage()
            }))) : t.itemsPerPage = i.itemsPerPage, n.$watch("totalItems", (function(e, r) {
              (angular.isDefined(e) || e !== r) && (n.totalPages = t.calculateTotalPages(), t.updatePage())
            }))
          }, t.calculateTotalPages = function() {
            var e = t.itemsPerPage < 1 ? 1 : Math.ceil(n.totalItems / t.itemsPerPage);
            return Math.max(e || 0, 1)
          }, t.render = function() {
            n.page = parseInt(t.ngModelCtrl.$viewValue, 10) || 1
          }, n.selectPage = function(e, r) {
            r && r.preventDefault(), (!n.ngDisabled || !r) && n.page !== e && e > 0 && e <= n.totalPages && (r && r.target && r.target.blur(), t
              .ngModelCtrl.$setViewValue(e), t.ngModelCtrl.$render())
          }, n.getText = function(e) {
            return n[e + "Text"] || t.config[e + "Text"]
          }, n.noPrevious = function() {
            return 1 === n.page
          }, n.noNext = function() {
            return n.page === n.totalPages
          }, t.updatePage = function() {
            t.setNumPages(n.$parent, n.totalPages), n.page > n.totalPages ? n.selectPage(n.totalPages) : t.ngModelCtrl.$render()
          }, n.$on("$destroy", (function() {
            for (; t._watchers.length;) t._watchers.shift()()
          }))
        }
      }
    }]), angular.module("ui.bootstrap.pager", ["ui.bootstrap.paging", "ui.bootstrap.tabindex"]).controller("UibPagerController", ["$scope", "$attrs",
      "uibPaging", "uibPagerConfig",
      function(e, t, n, r) {
        e.align = angular.isDefined(t.align) ? e.$parent.$eval(t.align) : r.align, n.create(this, e, t)
      }
    ]).constant("uibPagerConfig", {
      itemsPerPage: 10,
      previousText: "« Previous",
      nextText: "Next »",
      align: !0
    }).directive("uibPager", ["uibPagerConfig", function(e) {
      return {
        scope: {
          totalItems: "=",
          previousText: "@",
          nextText: "@",
          ngDisabled: "="
        },
        require: ["uibPager", "?ngModel"],
        restrict: "A",
        controller: "UibPagerController",
        controllerAs: "pager",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/pager/pager.html"
        },
        link: function(t, n, r, i) {
          n.addClass("pager");
          var o = i[0],
            a = i[1];
          a && o.init(a, e)
        }
      }
    }]), angular.module("ui.bootstrap.pagination", ["ui.bootstrap.paging", "ui.bootstrap.tabindex"]).controller("UibPaginationController", ["$scope",
      "$attrs", "$parse", "uibPaging", "uibPaginationConfig",
      function(e, t, n, r, i) {
        var o = this,
          a = angular.isDefined(t.maxSize) ? e.$parent.$eval(t.maxSize) : i.maxSize,
          s = angular.isDefined(t.rotate) ? e.$parent.$eval(t.rotate) : i.rotate,
          u = angular.isDefined(t.forceEllipses) ? e.$parent.$eval(t.forceEllipses) : i.forceEllipses,
          l = angular.isDefined(t.boundaryLinkNumbers) ? e.$parent.$eval(t.boundaryLinkNumbers) : i.boundaryLinkNumbers,
          c = angular.isDefined(t.pageLabel) ? function(n) {
            return e.$parent.$eval(t.pageLabel, {
              $page: n
            })
          } : angular.identity;

        function d(e, t, n) {
          return {
            number: e,
            text: t,
            active: n
          }
        }
        e.boundaryLinks = angular.isDefined(t.boundaryLinks) ? e.$parent.$eval(t.boundaryLinks) : i.boundaryLinks, e.directionLinks = angular
          .isDefined(t.directionLinks) ? e.$parent.$eval(t.directionLinks) : i.directionLinks, t.$set("role", "menu"), r.create(this, e, t), t
          .maxSize && o._watchers.push(e.$parent.$watch(n(t.maxSize), (function(e) {
            a = parseInt(e, 10), o.render()
          })));
        var p = this.render;
        this.render = function() {
          p(), e.page > 0 && e.page <= e.totalPages && (e.pages = function(e, t) {
            var n = [],
              r = 1,
              i = t,
              o = angular.isDefined(a) && a < t;
            o && (s ? (i = (r = Math.max(e - Math.floor(a / 2), 1)) + a - 1) > t && (r = (i = t) - a + 1) : (r = (Math.ceil(e / a) - 1) * a + 1,
              i = Math.min(r + a - 1, t)));
            for (var p = r; p <= i; p++) {
              var f = d(p, c(p), p === e);
              n.push(f)
            }
            if (o && a > 0 && (!s || u || l)) {
              if (r > 1) {
                if (!l || r > 3) {
                  var h = d(r - 1, "...", !1);
                  n.unshift(h)
                }
                if (l) {
                  if (3 === r) {
                    var g = d(2, "2", !1);
                    n.unshift(g)
                  }
                  var m = d(1, "1", !1);
                  n.unshift(m)
                }
              }
              if (i < t) {
                if (!l || i < t - 2) {
                  var v = d(i + 1, "...", !1);
                  n.push(v)
                }
                if (l) {
                  if (i === t - 2) {
                    var $ = d(t - 1, t - 1, !1);
                    n.push($)
                  }
                  var b = d(t, t, !1);
                  n.push(b)
                }
              }
            }
            return n
          }(e.page, e.totalPages))
        }
      }
    ]).constant("uibPaginationConfig", {
      itemsPerPage: 10,
      boundaryLinks: !1,
      boundaryLinkNumbers: !1,
      directionLinks: !0,
      firstText: "First",
      previousText: "Previous",
      nextText: "Next",
      lastText: "Last",
      rotate: !0,
      forceEllipses: !1
    }).directive("uibPagination", ["$parse", "uibPaginationConfig", function(e, t) {
      return {
        scope: {
          totalItems: "=",
          firstText: "@",
          previousText: "@",
          nextText: "@",
          lastText: "@",
          ngDisabled: "="
        },
        require: ["uibPagination", "?ngModel"],
        restrict: "A",
        controller: "UibPaginationController",
        controllerAs: "pagination",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/pagination/pagination.html"
        },
        link: function(e, n, r, i) {
          n.addClass("pagination");
          var o = i[0],
            a = i[1];
          a && o.init(a, t)
        }
      }
    }]), angular.module("ui.bootstrap.tooltip", ["ui.bootstrap.position", "ui.bootstrap.stackedMap"]).provider("$uibTooltip", (function() {
      var e = {
          placement: "top",
          placementClassPrefix: "",
          animation: !0,
          popupDelay: 0,
          popupCloseDelay: 0,
          useContentExp: !1
        },
        t = {
          mouseenter: "mouseleave",
          click: "click",
          outsideClick: "outsideClick",
          focus: "blur",
          none: ""
        },
        n = {};
      this.options = function(e) {
        angular.extend(n, e)
      }, this.setTriggers = function(e) {
        angular.extend(t, e)
      }, this.$get = ["$window", "$compile", "$timeout", "$document", "$uibPosition", "$interpolate", "$rootScope", "$parse", "$$stackedMap",
        function(r, i, o, a, s, u, l, c, d) {
          var p = d.createNew();

          function f(e) {
            if (27 === e.which) {
              var t = p.top();
              t && (t.value.close(), t = null)
            }
          }
          return a.on("keyup", f), l.$on("$destroy", (function() {
              a.off("keyup", f)
            })),
            function(r, l, d, f) {
              function h(e) {
                var n = (e || f.trigger || d).split(" "),
                  r = n.map((function(e) {
                    return t[e] || e
                  }));
                return {
                  show: n,
                  hide: r
                }
              }
              f = angular.extend({}, e, n, f);
              var g = r.replace(/[A-Z]/g, (function(e, t) {
                  return (t ? "-" : "") + e.toLowerCase()
                })),
                m = u.startSymbol(),
                v = u.endSymbol(),
                $ = "<div " + g + '-popup uib-title="' + m + "title" + v + '" ' + (f.useContentExp ? 'content-exp="contentExp()" ' : 'content="' +
                  m + "content" + v + '" ') + 'origin-scope="origScope" class="uib-position-measure ' + l +
                '" tooltip-animation-class="fade"uib-tooltip-classes ng-class="{ in: isOpen }" ></div>';
              return {
                compile: function(e, t) {
                  var n = i($);
                  return function(e, t, i, u) {
                    var d, g, m, v, $, b, y, w, x = !!angular.isDefined(f.appendToBody) && f.appendToBody,
                      C = h(void 0),
                      S = angular.isDefined(i[l + "Enable"]),
                      A = e.$new(!0),
                      T = !1,
                      k = !!angular.isDefined(i[l + "IsOpen"]) && c(i[l + "IsOpen"]),
                      D = !!f.useContentExp && c(i[r]),
                      M = [],
                      E = function() {
                        d && d.html() && (b || (b = o((function() {
                          var e = s.positionElements(t, d, A.placement, x),
                            n = angular.isDefined(d.offsetHeight) ? d.offsetHeight : d.prop("offsetHeight"),
                            r = x ? s.offset(t) : s.position(t);
                          d.css({
                            top: e.top + "px",
                            left: e.left + "px"
                          });
                          var i = e.placement.split("-");
                          d.hasClass(i[0]) || (d.removeClass(w.split("-")[0]), d.addClass(i[0])), d.hasClass(f
                              .placementClassPrefix + e.placement) || (d.removeClass(f.placementClassPrefix + w), d.addClass(f
                              .placementClassPrefix + e.placement)), y = o((function() {
                              var e = angular.isDefined(d.offsetHeight) ? d.offsetHeight : d.prop("offsetHeight"),
                                t = s.adjustTop(i, r, n, e);
                              t && d.css(t), y = null
                            }), 0, !1), d.hasClass("uib-position-measure") ? (s.positionArrow(d, e.placement), d.removeClass(
                              "uib-position-measure")) : w !== e.placement && s.positionArrow(d, e.placement), w = e.placement, b =
                            null
                        }), 0, !1)))
                      };

                    function O() {
                      A.isOpen ? P() : N()
                    }

                    function N() {
                      S && !e.$eval(i[l + "Enable"]) || (R(), function() {
                        A.title = i[l + "Title"], A.content = D ? D(e) : i[r];
                        A.popupClass = i[l + "Class"], A.placement = angular.isDefined(i[l + "Placement"]) ? i[l + "Placement"] : f
                          .placement;
                        var t = s.parsePlacement(A.placement);
                        w = t[1] ? t[0] + "-" + t[1] : t[0];
                        var n = parseInt(i[l + "PopupDelay"], 10),
                          o = parseInt(i[l + "PopupCloseDelay"], 10);
                        A.popupDelay = isNaN(n) ? f.popupDelay : n, A.popupCloseDelay = isNaN(o) ? f.popupCloseDelay : o
                      }(), A.popupDelay ? v || (v = o(q, A.popupDelay, !1)) : q())
                    }

                    function P() {
                      I(), A.popupCloseDelay ? $ || ($ = o(L, A.popupCloseDelay, !1)) : L()
                    }

                    function q() {
                      if (I(), R(), !A.content) return angular.noop;
                      ! function() {
                        if (d) return;
                        g = A.$new(), d = n(g, (function(e) {
                            x ? a.find("body").append(e) : t.after(e)
                          })), p.add(A, {
                            close: L
                          }),
                          function() {
                            M.length = 0, D ? (M.push(e.$watch(D, (function(e) {
                              A.content = e, !e && A.isOpen && L()
                            }))), M.push(g.$watch((function() {
                              T || (T = !0, g.$$postDigest((function() {
                                T = !1, A && A.isOpen && E()
                              })))
                            })))) : M.push(i.$observe(r, (function(e) {
                              A.content = e, !e && A.isOpen ? L() : E()
                            })));
                            M.push(i.$observe(l + "Title", (function(e) {
                              A.title = e, A.isOpen && E()
                            }))), M.push(i.$observe(l + "Placement", (function(e) {
                              A.placement = e || f.placement, A.isOpen && E()
                            })))
                          }()
                      }(), A.$evalAsync((function() {
                        A.isOpen = !0, _(!0), E()
                      }))
                    }

                    function I() {
                      v && (o.cancel(v), v = null), b && (o.cancel(b), b = null)
                    }

                    function L() {
                      A && A.$evalAsync((function() {
                        A && (A.isOpen = !1, _(!1), A.animation ? m || (m = o(j, 150, !1)) : j())
                      }))
                    }

                    function R() {
                      $ && (o.cancel($), $ = null), m && (o.cancel(m), m = null)
                    }

                    function j() {
                      I(), R(), M.length && (angular.forEach(M, (function(e) {
                        e()
                      })), M.length = 0), d && (d.remove(), d = null, y && o.cancel(y)), p.remove(A), g && (g.$destroy(), g = null)
                    }

                    function _(t) {
                      k && angular.isFunction(k.assign) && k.assign(e, t)
                    }

                    function V(e) {
                      A && A.isOpen && d && (t[0].contains(e.target) || d[0].contains(e.target) || P())
                    }

                    function U(e) {
                      27 === e.which && P()
                    }
                    A.origScope = e, A.isOpen = !1, A.contentExp = function() {
                      return A.content
                    }, i.$observe("disabled", (function(e) {
                      e && I(), e && A.isOpen && L()
                    })), k && e.$watch(k, (function(e) {
                      A && !e === A.isOpen && O()
                    }));
                    var H, F, B, G = function() {
                      C.show.forEach((function(e) {
                        "outsideClick" === e ? t.off("click", O) : (t.off(e, N), t.off(e, O)), t.off("keypress", U)
                      })), C.hide.forEach((function(e) {
                        "outsideClick" === e ? a.off("click", V) : t.off(e, P)
                      }))
                    };
                    H = [], F = [], B = e.$eval(i[l + "Trigger"]), G(), angular.isObject(B) ? (Object.keys(B).forEach((function(e) {
                      H.push(e), F.push(B[e])
                    })), C = {
                      show: H,
                      hide: F
                    }) : C = h(B), "none" !== C.show && C.show.forEach((function(e, n) {
                      "outsideClick" === e ? (t.on("click", O), a.on("click", V)) : e === C.hide[n] ? t.on(e, O) : e && (t.on(e, N), t
                        .on(C.hide[n], P)), t.on("keypress", U)
                    }));
                    var W, z = e.$eval(i[l + "Animation"]);
                    A.animation = angular.isDefined(z) ? !!z : f.animation;
                    var Y = l + "AppendToBody";
                    W = Y in i && void 0 === i[Y] || e.$eval(i[Y]), x = angular.isDefined(W) ? W : x, e.$on("$destroy", (function() {
                      G(), j(), A = null
                    }))
                  }
                }
              }
            }
        }
      ]
    })).directive("uibTooltipTemplateTransclude", ["$animate", "$sce", "$compile", "$templateRequest", function(e, t, n, r) {
      return {
        link: function(i, o, a) {
          var s, u, l, c = i.$eval(a.tooltipTemplateTranscludeScope),
            d = 0,
            p = function() {
              u && (u.remove(), u = null), s && (s.$destroy(), s = null), l && (e.leave(l).then((function() {
                u = null
              })), u = l, l = null)
            };
          i.$watch(t.parseAsResourceUrl(a.uibTooltipTemplateTransclude), (function(t) {
            var a = ++d;
            t ? (r(t, !0).then((function(r) {
              if (a === d) {
                var i = c.$new(),
                  u = n(r)(i, (function(t) {
                    p(), e.enter(t, o)
                  }));
                l = u, (s = i).$emit("$includeContentLoaded", t)
              }
            }), (function() {
              a === d && (p(), i.$emit("$includeContentError", t))
            })), i.$emit("$includeContentRequested", t)) : p()
          })), i.$on("$destroy", p)
        }
      }
    }]).directive("uibTooltipClasses", ["$uibPosition", function(e) {
      return {
        restrict: "A",
        link: function(t, n, r) {
          if (t.placement) {
            var i = e.parsePlacement(t.placement);
            n.addClass(i[0])
          }
          t.popupClass && n.addClass(t.popupClass), t.animation && n.addClass(r.tooltipAnimationClass)
        }
      }
    }]).directive("uibTooltipPopup", (function() {
      return {
        restrict: "A",
        scope: {
          content: "@"
        },
        templateUrl: "uib/template/tooltip/tooltip-popup.html"
      }
    })).directive("uibTooltip", ["$uibTooltip", function(e) {
      return e("uibTooltip", "tooltip", "mouseenter")
    }]).directive("uibTooltipTemplatePopup", (function() {
      return {
        restrict: "A",
        scope: {
          contentExp: "&",
          originScope: "&"
        },
        templateUrl: "uib/template/tooltip/tooltip-template-popup.html"
      }
    })).directive("uibTooltipTemplate", ["$uibTooltip", function(e) {
      return e("uibTooltipTemplate", "tooltip", "mouseenter", {
        useContentExp: !0
      })
    }]).directive("uibTooltipHtmlPopup", (function() {
      return {
        restrict: "A",
        scope: {
          contentExp: "&"
        },
        templateUrl: "uib/template/tooltip/tooltip-html-popup.html"
      }
    })).directive("uibTooltipHtml", ["$uibTooltip", function(e) {
      return e("uibTooltipHtml", "tooltip", "mouseenter", {
        useContentExp: !0
      })
    }]), angular.module("ui.bootstrap.popover", ["ui.bootstrap.tooltip"]).directive("uibPopoverTemplatePopup", (function() {
      return {
        restrict: "A",
        scope: {
          uibTitle: "@",
          contentExp: "&",
          originScope: "&"
        },
        templateUrl: "uib/template/popover/popover-template.html"
      }
    })).directive("uibPopoverTemplate", ["$uibTooltip", function(e) {
      return e("uibPopoverTemplate", "popover", "click", {
        useContentExp: !0
      })
    }]).directive("uibPopoverHtmlPopup", (function() {
      return {
        restrict: "A",
        scope: {
          contentExp: "&",
          uibTitle: "@"
        },
        templateUrl: "uib/template/popover/popover-html.html"
      }
    })).directive("uibPopoverHtml", ["$uibTooltip", function(e) {
      return e("uibPopoverHtml", "popover", "click", {
        useContentExp: !0
      })
    }]).directive("uibPopoverPopup", (function() {
      return {
        restrict: "A",
        scope: {
          uibTitle: "@",
          content: "@"
        },
        templateUrl: "uib/template/popover/popover.html"
      }
    })).directive("uibPopover", ["$uibTooltip", function(e) {
      return e("uibPopover", "popover", "click")
    }]), angular.module("ui.bootstrap.progressbar", []).constant("uibProgressConfig", {
      animate: !0,
      max: 100
    }).controller("UibProgressController", ["$scope", "$attrs", "uibProgressConfig", function(e, t, n) {
      var r = this,
        i = angular.isDefined(t.animate) ? e.$parent.$eval(t.animate) : n.animate;

      function o() {
        return angular.isDefined(e.maxParam) ? e.maxParam : n.max
      }
      this.bars = [], e.max = o(), this.addBar = function(e, t, n) {
        i || t.css({
          transition: "none"
        }), this.bars.push(e), e.max = o(), e.title = n && angular.isDefined(n.title) ? n.title : "progressbar", e.$watch("value", (function(
        t) {
          e.recalculatePercentage()
        })), e.recalculatePercentage = function() {
          var t = r.bars.reduce((function(e, t) {
            return t.percent = +(100 * t.value / t.max).toFixed(2), e + t.percent
          }), 0);
          t > 100 && (e.percent -= t - 100)
        }, e.$on("$destroy", (function() {
          t = null, r.removeBar(e)
        }))
      }, this.removeBar = function(e) {
        this.bars.splice(this.bars.indexOf(e), 1), this.bars.forEach((function(e) {
          e.recalculatePercentage()
        }))
      }, e.$watch("maxParam", (function(e) {
        r.bars.forEach((function(e) {
          e.max = o(), e.recalculatePercentage()
        }))
      }))
    }]).directive("uibProgress", (function() {
      return {
        replace: !0,
        transclude: !0,
        controller: "UibProgressController",
        require: "uibProgress",
        scope: {
          maxParam: "=?max"
        },
        templateUrl: "uib/template/progressbar/progress.html"
      }
    })).directive("uibBar", (function() {
      return {
        replace: !0,
        transclude: !0,
        require: "^uibProgress",
        scope: {
          value: "=",
          type: "@"
        },
        templateUrl: "uib/template/progressbar/bar.html",
        link: function(e, t, n, r) {
          r.addBar(e, t, n)
        }
      }
    })).directive("uibProgressbar", (function() {
      return {
        replace: !0,
        transclude: !0,
        controller: "UibProgressController",
        scope: {
          value: "=",
          maxParam: "=?max",
          type: "@"
        },
        templateUrl: "uib/template/progressbar/progressbar.html",
        link: function(e, t, n, r) {
          r.addBar(e, angular.element(t.children()[0]), {
            title: n.title
          })
        }
      }
    })), angular.module("ui.bootstrap.rating", []).constant("uibRatingConfig", {
      max: 5,
      stateOn: null,
      stateOff: null,
      enableReset: !0,
      titles: ["one", "two", "three", "four", "five"]
    }).controller("UibRatingController", ["$scope", "$attrs", "uibRatingConfig", function(e, t, n) {
      var r = {
          $setViewValue: angular.noop
        },
        i = this;
      this.init = function(i) {
        (r = i).$render = this.render, r.$formatters.push((function(e) {
            return angular.isNumber(e) && e << 0 !== e && (e = Math.round(e)), e
          })), this.stateOn = angular.isDefined(t.stateOn) ? e.$parent.$eval(t.stateOn) : n.stateOn, this.stateOff = angular.isDefined(t
          .stateOff) ? e.$parent.$eval(t.stateOff) : n.stateOff, this.enableReset = angular.isDefined(t.enableReset) ? e.$parent.$eval(t
            .enableReset) : n.enableReset;
        var o = angular.isDefined(t.titles) ? e.$parent.$eval(t.titles) : n.titles;
        this.titles = angular.isArray(o) && o.length > 0 ? o : n.titles;
        var a = angular.isDefined(t.ratingStates) ? e.$parent.$eval(t.ratingStates) : new Array(angular.isDefined(t.max) ? e.$parent.$eval(t
          .max) : n.max);
        e.range = this.buildTemplateObjects(a)
      }, this.buildTemplateObjects = function(e) {
        for (var t = 0, n = e.length; t < n; t++) e[t] = angular.extend({
          index: t
        }, {
          stateOn: this.stateOn,
          stateOff: this.stateOff,
          title: this.getTitle(t)
        }, e[t]);
        return e
      }, this.getTitle = function(e) {
        return e >= this.titles.length ? e + 1 : this.titles[e]
      }, e.rate = function(t) {
        if (!e.readonly && t >= 0 && t <= e.range.length) {
          var n = i.enableReset && r.$viewValue === t ? 0 : t;
          r.$setViewValue(n), r.$render()
        }
      }, e.enter = function(t) {
        e.readonly || (e.value = t), e.onHover({
          value: t
        })
      }, e.reset = function() {
        e.value = r.$viewValue, e.onLeave()
      }, e.onKeydown = function(t) {
        /(37|38|39|40)/.test(t.which) && (t.preventDefault(), t.stopPropagation(), e.rate(e.value + (38 === t.which || 39 === t.which ? 1 : -1)))
      }, this.render = function() {
        e.value = r.$viewValue, e.title = i.getTitle(e.value - 1)
      }
    }]).directive("uibRating", (function() {
      return {
        require: ["uibRating", "ngModel"],
        restrict: "A",
        scope: {
          readonly: "=?readOnly",
          onHover: "&",
          onLeave: "&"
        },
        controller: "UibRatingController",
        templateUrl: "uib/template/rating/rating.html",
        link: function(e, t, n, r) {
          var i = r[0],
            o = r[1];
          i.init(o)
        }
      }
    })), angular.module("ui.bootstrap.tabs", []).controller("UibTabsetController", ["$scope", function(e) {
      var t, n, r = this;

      function i(e) {
        for (var t = 0; t < r.tabs.length; t++)
          if (r.tabs[t].index === e) return t
      }
      r.tabs = [], r.select = function(e, o) {
        if (!n) {
          var a = i(t),
            s = r.tabs[a];
          if (s) {
            if (s.tab.onDeselect({
                $event: o,
                $selectedIndex: e
              }), o && o.isDefaultPrevented()) return;
            s.tab.active = !1
          }
          var u = r.tabs[e];
          u ? (u.tab.onSelect({
            $event: o
          }), u.tab.active = !0, r.active = u.index, t = u.index) : !u && angular.isDefined(t) && (r.active = null, t = null)
        }
      }, r.addTab = function(e) {
        if (r.tabs.push({
            tab: e,
            index: e.index
          }), r.tabs.sort((function(e, t) {
            return e.index > t.index ? 1 : e.index < t.index ? -1 : 0
          })), e.index === r.active || !angular.isDefined(r.active) && 1 === r.tabs.length) {
          var t = i(e.index);
          r.select(t)
        }
      }, r.removeTab = function(e) {
        for (var t, n = 0; n < r.tabs.length; n++)
          if (r.tabs[n].tab === e) {
            t = n;
            break
          } if (r.tabs[t].index === r.active) {
          var i = t === r.tabs.length - 1 ? t - 1 : t + 1 % r.tabs.length;
          r.select(i)
        }
        r.tabs.splice(t, 1)
      }, e.$watch("tabset.active", (function(e) {
        angular.isDefined(e) && e !== t && r.select(i(e))
      })), e.$on("$destroy", (function() {
        n = !0
      }))
    }]).directive("uibTabset", (function() {
      return {
        transclude: !0,
        replace: !0,
        scope: {},
        bindToController: {
          active: "=?",
          type: "@"
        },
        controller: "UibTabsetController",
        controllerAs: "tabset",
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/tabs/tabset.html"
        },
        link: function(e, t, n) {
          e.vertical = !!angular.isDefined(n.vertical) && e.$parent.$eval(n.vertical), e.justified = !!angular.isDefined(n.justified) && e.$parent
            .$eval(n.justified)
        }
      }
    })).directive("uibTab", ["$parse", function(e) {
      return {
        require: "^uibTabset",
        replace: !0,
        templateUrl: function(e, t) {
          return t.templateUrl || "uib/template/tabs/tab.html"
        },
        transclude: !0,
        scope: {
          heading: "@",
          index: "=?",
          classes: "@?",
          onSelect: "&select",
          onDeselect: "&deselect"
        },
        controller: function() {},
        controllerAs: "tab",
        link: function(t, n, r, i, o) {
          t.disabled = !1, r.disable && t.$parent.$watch(e(r.disable), (function(e) {
            t.disabled = !!e
          })), angular.isUndefined(r.index) && (i.tabs && i.tabs.length ? t.index = Math.max.apply(null, i.tabs.map((function(e) {
            return e.index
          }))) + 1 : t.index = 0), angular.isUndefined(r.classes) && (t.classes = ""), t.select = function(e) {
            if (!t.disabled) {
              for (var n, r = 0; r < i.tabs.length; r++)
                if (i.tabs[r].tab === t) {
                  n = r;
                  break
                } i.select(n, e)
            }
          }, i.addTab(t), t.$on("$destroy", (function() {
            i.removeTab(t)
          })), t.$transcludeFn = o
        }
      }
    }]).directive("uibTabHeadingTransclude", (function() {
      return {
        restrict: "A",
        require: "^uibTab",
        link: function(e, t) {
          e.$watch("headingElement", (function(e) {
            e && (t.html(""), t.append(e))
          }))
        }
      }
    })).directive("uibTabContentTransclude", (function() {
      return {
        restrict: "A",
        require: "^uibTabset",
        link: function(e, t, n) {
          var r = e.$eval(n.uibTabContentTransclude).tab;
          r.$transcludeFn(r.$parent, (function(e) {
            angular.forEach(e, (function(e) {
              ! function(e) {
                return e.tagName && (e.hasAttribute("uib-tab-heading") || e.hasAttribute("data-uib-tab-heading") || e.hasAttribute(
                    "x-uib-tab-heading") || "uib-tab-heading" === e.tagName.toLowerCase() || "data-uib-tab-heading" === e.tagName
                  .toLowerCase() || "x-uib-tab-heading" === e.tagName.toLowerCase() || "uib:tab-heading" === e.tagName.toLowerCase()
                  )
              }(e) ? t.append(e): r.headingElement = e
            }))
          }))
        }
      }
    })), angular.module("ui.bootstrap.timepicker", []).constant("uibTimepickerConfig", {
      hourStep: 1,
      minuteStep: 1,
      secondStep: 1,
      showMeridian: !0,
      showSeconds: !1,
      meridians: null,
      readonlyInput: !1,
      mousewheel: !0,
      arrowkeys: !0,
      showSpinners: !0,
      templateUrl: "uib/template/timepicker/timepicker.html"
    }).controller("UibTimepickerController", ["$scope", "$element", "$attrs", "$parse", "$log", "$locale", "uibTimepickerConfig", function(e, t, n, r,
      i, o, a) {
      var s, u, l, c = new Date,
        d = [],
        p = {
          $setViewValue: angular.noop
        },
        f = angular.isDefined(n.meridians) ? e.$parent.$eval(n.meridians) : a.meridians || o.DATETIME_FORMATS.AMPMS,
        h = !angular.isDefined(n.padHours) || e.$parent.$eval(n.padHours);
      e.tabindex = angular.isDefined(n.tabindex) ? n.tabindex : 0, t.removeAttr("tabindex"), this.init = function(t, r) {
        (p = t).$render = this.render, p.$formatters.unshift((function(e) {
          return e ? new Date(e) : null
        }));
        var i = r.eq(0),
          o = r.eq(1),
          c = r.eq(2);
        s = i.controller("ngModel"), u = o.controller("ngModel"), l = c.controller("ngModel"), (angular.isDefined(n.mousewheel) ? e.$parent.$eval(
          n.mousewheel) : a.mousewheel) && this.setupMousewheelEvents(i, o, c), (angular.isDefined(n.arrowkeys) ? e.$parent.$eval(n.arrowkeys) :
          a.arrowkeys) && this.setupArrowkeyEvents(i, o, c), e.readonlyInput = angular.isDefined(n.readonlyInput) ? e.$parent.$eval(n
          .readonlyInput) : a.readonlyInput, this.setupInputEvents(i, o, c)
      };
      var g = a.hourStep;
      n.hourStep && d.push(e.$parent.$watch(r(n.hourStep), (function(e) {
        g = +e
      })));
      var m, v, $ = a.minuteStep;
      n.minuteStep && d.push(e.$parent.$watch(r(n.minuteStep), (function(e) {
        $ = +e
      }))), d.push(e.$parent.$watch(r(n.min), (function(e) {
        var t = new Date(e);
        m = isNaN(t) ? void 0 : t
      }))), d.push(e.$parent.$watch(r(n.max), (function(e) {
        var t = new Date(e);
        v = isNaN(t) ? void 0 : t
      })));
      var b = !1;
      n.ngDisabled && d.push(e.$parent.$watch(r(n.ngDisabled), (function(e) {
        b = e
      }))), e.noIncrementHours = function() {
        var e = D(c, 60 * g);
        return b || e > v || e < c && e < m
      }, e.noDecrementHours = function() {
        var e = D(c, 60 * -g);
        return b || e < m || e > c && e > v
      }, e.noIncrementMinutes = function() {
        var e = D(c, $);
        return b || e > v || e < c && e < m
      }, e.noDecrementMinutes = function() {
        var e = D(c, -$);
        return b || e < m || e > c && e > v
      }, e.noIncrementSeconds = function() {
        var e = M(c, y);
        return b || e > v || e < c && e < m
      }, e.noDecrementSeconds = function() {
        var e = M(c, -y);
        return b || e < m || e > c && e > v
      }, e.noToggleMeridian = function() {
        return c.getHours() < 12 ? b || D(c, 720) > v : b || D(c, -720) < m
      };
      var y = a.secondStep;

      function w() {
        var t = +e.hours;
        if ((e.showMeridian ? t > 0 && t < 13 : t >= 0 && t < 24) && "" !== e.hours) return e.showMeridian && (12 === t && (t = 0), e.meridian ===
          f[1] && (t += 12)), t
      }

      function x() {
        var t = +e.minutes;
        if (t >= 0 && t < 60 && "" !== e.minutes) return t
      }

      function C(e, t) {
        return null === e ? "" : angular.isDefined(e) && e.toString().length < 2 && !t ? "0" + e : e.toString()
      }

      function S(e) {
        A(), p.$setViewValue(new Date(c)), T(e)
      }

      function A() {
        s && s.$setValidity("hours", !0), u && u.$setValidity("minutes", !0), l && l.$setValidity("seconds", !0), p.$setValidity("time", !0), e
          .invalidHours = !1, e.invalidMinutes = !1, e.invalidSeconds = !1
      }

      function T(t) {
        if (p.$modelValue) {
          var n = c.getHours(),
            r = c.getMinutes(),
            i = c.getSeconds();
          e.showMeridian && (n = 0 === n || 12 === n ? 12 : n % 12), e.hours = "h" === t ? n : C(n, !h), "m" !== t && (e.minutes = C(r)), e
            .meridian = c.getHours() < 12 ? f[0] : f[1], "s" !== t && (e.seconds = C(i)), e.meridian = c.getHours() < 12 ? f[0] : f[1]
        } else e.hours = null, e.minutes = null, e.seconds = null, e.meridian = f[0]
      }

      function k(e) {
        c = M(c, e), S()
      }

      function D(e, t) {
        return M(e, 60 * t)
      }

      function M(e, t) {
        var n = new Date(e.getTime() + 1e3 * t),
          r = new Date(e);
        return r.setHours(n.getHours(), n.getMinutes(), n.getSeconds()), r
      }

      function E() {
        return (null === e.hours || "" === e.hours) && (null === e.minutes || "" === e.minutes) && (!e.showSeconds || e.showSeconds && (null === e
          .seconds || "" === e.seconds))
      }
      n.secondStep && d.push(e.$parent.$watch(r(n.secondStep), (function(e) {
        y = +e
      }))), e.showSeconds = a.showSeconds, n.showSeconds && d.push(e.$parent.$watch(r(n.showSeconds), (function(t) {
        e.showSeconds = !!t
      }))), e.showMeridian = a.showMeridian, n.showMeridian && d.push(e.$parent.$watch(r(n.showMeridian), (function(t) {
        if (e.showMeridian = !!t, p.$error.time) {
          var n = w(),
            r = x();
          angular.isDefined(n) && angular.isDefined(r) && (c.setHours(n), S())
        } else T()
      }))), this.setupMousewheelEvents = function(t, n, r) {
        var i = function(e) {
          e.originalEvent && (e = e.originalEvent);
          var t = e.wheelDelta ? e.wheelDelta : -e.deltaY;
          return e.detail || t > 0
        };
        t.on("mousewheel wheel", (function(t) {
          b || e.$apply(i(t) ? e.incrementHours() : e.decrementHours()), t.preventDefault()
        })), n.on("mousewheel wheel", (function(t) {
          b || e.$apply(i(t) ? e.incrementMinutes() : e.decrementMinutes()), t.preventDefault()
        })), r.on("mousewheel wheel", (function(t) {
          b || e.$apply(i(t) ? e.incrementSeconds() : e.decrementSeconds()), t.preventDefault()
        }))
      }, this.setupArrowkeyEvents = function(t, n, r) {
        t.on("keydown", (function(t) {
          b || (38 === t.which ? (t.preventDefault(), e.incrementHours(), e.$apply()) : 40 === t.which && (t.preventDefault(), e
            .decrementHours(), e.$apply()))
        })), n.on("keydown", (function(t) {
          b || (38 === t.which ? (t.preventDefault(), e.incrementMinutes(), e.$apply()) : 40 === t.which && (t.preventDefault(), e
            .decrementMinutes(), e.$apply()))
        })), r.on("keydown", (function(t) {
          b || (38 === t.which ? (t.preventDefault(), e.incrementSeconds(), e.$apply()) : 40 === t.which && (t.preventDefault(), e
            .decrementSeconds(), e.$apply()))
        }))
      }, this.setupInputEvents = function(t, n, r) {
        if (e.readonlyInput) return e.updateHours = angular.noop, e.updateMinutes = angular.noop, void(e.updateSeconds = angular.noop);
        var i = function(t, n, r) {
          p.$setViewValue(null), p.$setValidity("time", !1), angular.isDefined(t) && (e.invalidHours = t, s && s.$setValidity("hours", !1)),
            angular.isDefined(n) && (e.invalidMinutes = n, u && u.$setValidity("minutes", !1)), angular.isDefined(r) && (e.invalidSeconds = r,
              l && l.$setValidity("seconds", !1))
        };
        e.updateHours = function() {
          var e = w(),
            t = x();
          p.$setDirty(), angular.isDefined(e) && angular.isDefined(t) ? (c.setHours(e), c.setMinutes(t), c < m || c > v ? i(!0) : S("h")) : i(!
            0)
        }, t.on("blur", (function(t) {
          p.$setTouched(), E() ? A() : null === e.hours || "" === e.hours ? i(!0) : !e.invalidHours && e.hours < 10 && e.$apply((
        function() {
            e.hours = C(e.hours, !h)
          }))
        })), e.updateMinutes = function() {
          var e = x(),
            t = w();
          p.$setDirty(), angular.isDefined(e) && angular.isDefined(t) ? (c.setHours(t), c.setMinutes(e), c < m || c > v ? i(void 0, !0) : S(
            "m")) : i(void 0, !0)
        }, n.on("blur", (function(t) {
          p.$setTouched(), E() ? A() : null === e.minutes ? i(void 0, !0) : !e.invalidMinutes && e.minutes < 10 && e.$apply((function() {
            e.minutes = C(e.minutes)
          }))
        })), e.updateSeconds = function() {
          var t = function() {
            var t = +e.seconds;
            return t >= 0 && t < 60 ? t : void 0
          }();
          p.$setDirty(), angular.isDefined(t) ? (c.setSeconds(t), S("s")) : i(void 0, void 0, !0)
        }, r.on("blur", (function(t) {
          E() ? A() : !e.invalidSeconds && e.seconds < 10 && e.$apply((function() {
            e.seconds = C(e.seconds)
          }))
        }))
      }, this.render = function() {
        var t = p.$viewValue;
        isNaN(t) ? (p.$setValidity("time", !1), i.error(
          'Timepicker directive: "ng-model" value must be a Date object, a number of milliseconds since 01.01.1970 or a string representing an RFC2822 or ISO 8601 date.'
          )) : (t && (c = t), c < m || c > v ? (p.$setValidity("time", !1), e.invalidHours = !0, e.invalidMinutes = !0) : A(), T())
      }, e.showSpinners = angular.isDefined(n.showSpinners) ? e.$parent.$eval(n.showSpinners) : a.showSpinners, e.incrementHours = function() {
        e.noIncrementHours() || k(60 * g * 60)
      }, e.decrementHours = function() {
        e.noDecrementHours() || k(60 * -g * 60)
      }, e.incrementMinutes = function() {
        e.noIncrementMinutes() || k(60 * $)
      }, e.decrementMinutes = function() {
        e.noDecrementMinutes() || k(60 * -$)
      }, e.incrementSeconds = function() {
        e.noIncrementSeconds() || k(y)
      }, e.decrementSeconds = function() {
        e.noDecrementSeconds() || k(-y)
      }, e.toggleMeridian = function() {
        var t = x(),
          n = w();
        e.noToggleMeridian() || (angular.isDefined(t) && angular.isDefined(n) ? k(720 * (c.getHours() < 12 ? 60 : -60)) : e.meridian = e
          .meridian === f[0] ? f[1] : f[0])
      }, e.blur = function() {
        p.$setTouched()
      }, e.$on("$destroy", (function() {
        for (; d.length;) d.shift()()
      }))
    }]).directive("uibTimepicker", ["uibTimepickerConfig", function(e) {
      return {
        require: ["uibTimepicker", "?^ngModel"],
        restrict: "A",
        controller: "UibTimepickerController",
        controllerAs: "timepicker",
        scope: {},
        templateUrl: function(t, n) {
          return n.templateUrl || e.templateUrl
        },
        link: function(e, t, n, r) {
          var i = r[0],
            o = r[1];
          o && i.init(o, t.find("input"))
        }
      }
    }]), angular.module("ui.bootstrap.typeahead", ["ui.bootstrap.debounce", "ui.bootstrap.position"]).factory("uibTypeaheadParser", ["$parse", function(
      e) {
      var t = /^\s*([\s\S]+?)(?:\s+as\s+([\s\S]+?))?\s+for\s+(?:([\$\w][\$\w\d]*))\s+in\s+([\s\S]+?)$/;
      return {
        parse: function(n) {
          var r = n.match(t);
          if (!r) throw new Error(
            'Expected typeahead specification in form of "_modelValue_ (as _label_)? for _item_ in _collection_" but got "' + n + '".');
          return {
            itemName: r[3],
            source: e(r[4]),
            viewMapper: e(r[2] || r[1]),
            modelMapper: e(r[1])
          }
        }
      }
    }]).controller("UibTypeaheadController", ["$scope", "$element", "$attrs", "$compile", "$parse", "$q", "$timeout", "$document", "$window",
      "$rootScope", "$$debounce", "$uibPosition", "uibTypeaheadParser",
      function(e, t, n, r, i, o, a, s, u, l, c, d, p) {
        var f, h, g = [9, 13, 27, 38, 40],
          m = e.$eval(n.typeaheadMinLength);
        m || 0 === m || (m = 1), e.$watch(n.typeaheadMinLength, (function(e) {
          m = e || 0 === e ? e : 1
        }));
        var v = e.$eval(n.typeaheadWaitMs) || 0,
          $ = !1 !== e.$eval(n.typeaheadEditable);
        e.$watch(n.typeaheadEditable, (function(e) {
          $ = !1 !== e
        }));
        var b, y, w = i(n.typeaheadLoading).assign || angular.noop,
          x = n.typeaheadShouldSelect ? i(n.typeaheadShouldSelect) : function(e, t) {
            var n = t.$event;
            return 13 === n.which || 9 === n.which
          },
          C = i(n.typeaheadOnSelect),
          S = !!angular.isDefined(n.typeaheadSelectOnBlur) && e.$eval(n.typeaheadSelectOnBlur),
          A = i(n.typeaheadNoResults).assign || angular.noop,
          T = n.typeaheadInputFormatter ? i(n.typeaheadInputFormatter) : void 0,
          k = !!n.typeaheadAppendToBody && e.$eval(n.typeaheadAppendToBody),
          D = n.typeaheadAppendTo ? e.$eval(n.typeaheadAppendTo) : null,
          M = !1 !== e.$eval(n.typeaheadFocusFirst),
          E = !!n.typeaheadSelectOnExact && e.$eval(n.typeaheadSelectOnExact),
          O = i(n.typeaheadIsOpen).assign || angular.noop,
          N = e.$eval(n.typeaheadShowHint) || !1,
          P = i(n.ngModel),
          q = i(n.ngModel + "($$$p)"),
          I = p.parse(n.uibTypeahead),
          L = e.$new(),
          R = e.$on("$destroy", (function() {
            L.$destroy()
          }));
        L.$on("$destroy", R);
        var j, _, V = "typeahead-" + L.$id + "-" + Math.floor(1e4 * Math.random());
        t.attr({
          "aria-autocomplete": "list",
          "aria-expanded": !1,
          "aria-owns": V
        }), N && ((j = angular.element("<div></div>")).css("position", "relative"), t.after(j), (_ = t.clone()).attr("placeholder", ""), _.attr(
          "tabindex", "-1"), _.val(""), _.css({
          position: "absolute",
          top: "0px",
          left: "0px",
          "border-color": "transparent",
          "box-shadow": "none",
          opacity: 1,
          background: "none 0% 0% / auto repeat scroll padding-box border-box rgb(255, 255, 255)",
          color: "#999"
        }), t.css({
          position: "relative",
          "vertical-align": "top",
          "background-color": "transparent"
        }), _.attr("id") && _.removeAttr("id"), j.append(_), _.after(t));
        var U = angular.element("<div uib-typeahead-popup></div>");
        U.attr({
          id: V,
          matches: "matches",
          active: "activeIdx",
          select: "select(activeIdx, evt)",
          "move-in-progress": "moveInProgress",
          query: "query",
          position: "position",
          "assign-is-open": "assignIsOpen(isOpen)",
          debounce: "debounceUpdate"
        }), angular.isDefined(n.typeaheadTemplateUrl) && U.attr("template-url", n.typeaheadTemplateUrl), angular.isDefined(n
          .typeaheadPopupTemplateUrl) && U.attr("popup-template-url", n.typeaheadPopupTemplateUrl);
        var H = function() {
            L.matches = [], L.activeIdx = -1, t.attr("aria-expanded", !1), N && _.val("")
          },
          F = function(e) {
            return V + "-option-" + e
          };
        L.$watch("activeIdx", (function(e) {
          e < 0 ? t.removeAttr("aria-activedescendant") : t.attr("aria-activedescendant", F(e))
        }));
        var B = function(n, r) {
          var i = {
            $viewValue: n
          };
          w(e, !0), A(e, !1), o.when(I.source(e, i)).then((function(o) {
            var a = n === f.$viewValue;
            if (a && b)
              if (o && o.length > 0) {
                L.activeIdx = M ? 0 : -1, A(e, !1), L.matches.length = 0;
                for (var s = 0; s < o.length; s++) i[I.itemName] = o[s], L.matches.push({
                  id: F(s),
                  label: I.viewMapper(L, i),
                  model: o[s]
                });
                if (L.query = n, Y(), t.attr("aria-expanded", !0), E && 1 === L.matches.length && function(e, t) {
                    return !!(L.matches.length > t && e) && e.toUpperCase() === L.matches[t].label.toUpperCase()
                  }(n, 0) && (angular.isNumber(L.debounceUpdate) || angular.isObject(L.debounceUpdate) ? c((function() {
                    L.select(0, r)
                  }), angular.isNumber(L.debounceUpdate) ? L.debounceUpdate : L.debounceUpdate.default) : L.select(0, r)), N) {
                  var u = L.matches[0].label;
                  angular.isString(n) && n.length > 0 && u.slice(0, n.length).toUpperCase() === n.toUpperCase() ? _.val(n + u.slice(n
                    .length)) : _.val("")
                }
              } else H(), A(e, !0);
            a && w(e, !1)
          }), (function() {
            H(), w(e, !1), A(e, !0)
          }))
        };
        k && (angular.element(u).on("resize", z), s.find("body").on("scroll", z));
        var G, W = c((function() {
          L.matches.length && Y(), L.moveInProgress = !1
        }), 200);

        function z() {
          L.moveInProgress || (L.moveInProgress = !0, L.$digest()), W()
        }

        function Y() {
          L.position = k ? d.offset(t) : d.position(t), L.position.top += t.prop("offsetHeight")
        }
        L.moveInProgress = !1, L.query = void 0;
        var K = function() {
          G && a.cancel(G)
        };
        H(), L.assignIsOpen = function(t) {
          O(e, t)
        }, L.select = function(r, i) {
          var o, s, u = {};
          y = !0, u[I.itemName] = s = L.matches[r].model, o = I.modelMapper(e, u),
            function(t, n) {
              angular.isFunction(P(e)) && h.getOption("getterSetter") ? q(t, {
                $$$p: n
              }) : P.assign(t, n)
            }(e, o), f.$setValidity("editable", !0), f.$setValidity("parse", !0), C(e, {
              $item: s,
              $model: o,
              $label: I.viewMapper(e, u),
              $event: i
            }), H(), !1 !== L.$eval(n.typeaheadFocusOnSelect) && a((function() {
              t[0].focus()
            }), 0, !1)
        }, t.on("keydown", (function(t) {
          if (0 !== L.matches.length && -1 !== g.indexOf(t.which)) {
            var n, r = x(e, {
              $event: t
            });
            if (-1 === L.activeIdx && r || 9 === t.which && t.shiftKey) return H(), void L.$digest();
            switch (t.preventDefault(), t.which) {
              case 27:
                t.stopPropagation(), H(), e.$digest();
                break;
              case 38:
                L.activeIdx = (L.activeIdx > 0 ? L.activeIdx : L.matches.length) - 1, L.$digest(), (n = U[0].querySelectorAll(
                  ".uib-typeahead-match")[L.activeIdx]).parentNode.scrollTop = n.offsetTop;
                break;
              case 40:
                L.activeIdx = (L.activeIdx + 1) % L.matches.length, L.$digest(), (n = U[0].querySelectorAll(".uib-typeahead-match")[L
                  .activeIdx]).parentNode.scrollTop = n.offsetTop;
                break;
              default:
                r && L.$apply((function() {
                  angular.isNumber(L.debounceUpdate) || angular.isObject(L.debounceUpdate) ? c((function() {
                    L.select(L.activeIdx, t)
                  }), angular.isNumber(L.debounceUpdate) ? L.debounceUpdate : L.debounceUpdate.default) : L.select(L.activeIdx, t)
                }))
            }
          }
        })), t.on("focus", (function(e) {
          b = !0, 0 !== m || f.$viewValue || a((function() {
            B(f.$viewValue, e)
          }), 0)
        })), t.on("blur", (function(e) {
          S && L.matches.length && -1 !== L.activeIdx && !y && (y = !0, L.$apply((function() {
            angular.isObject(L.debounceUpdate) && angular.isNumber(L.debounceUpdate.blur) ? c((function() {
              L.select(L.activeIdx, e)
            }), L.debounceUpdate.blur) : L.select(L.activeIdx, e)
          }))), !$ && f.$error.editable && (f.$setViewValue(), L.$apply((function() {
            f.$setValidity("editable", !0), f.$setValidity("parse", !0)
          })), t.val("")), b = !1, y = !1
        }));
        var Q = function(n) {
          t[0] !== n.target && 3 !== n.which && 0 !== L.matches.length && (H(), l.$$phase || e.$digest())
        };
        s.on("click", Q), e.$on("$destroy", (function() {
          s.off("click", Q), (k || D) && X.remove(), k && (angular.element(u).off("resize", z), s.find("body").off("scroll", z)), U.remove(),
            N && j.remove()
        }));
        var X = r(U)(L);
        k ? s.find("body").append(X) : D ? angular.element(D).eq(0).append(X) : t.after(X), this.init = function(t) {
          h = function(e) {
            var t;
            angular.version.minor < 6 ? (t = e.$options || {}).getOption = function(e) {
              return t[e]
            } : t = e.$options;
            return t
          }(f = t), L.debounceUpdate = i(h.getOption("debounce"))(e), f.$parsers.unshift((function(t) {
            return b = !0, 0 === m || t && t.length >= m ? v > 0 ? (K(), function(e) {
              G = a((function() {
                B(e)
              }), v)
            }(t)) : B(t) : (w(e, !1), K(), H()), $ ? t : t ? void f.$setValidity("editable", !1) : (f.$setValidity("editable", !0), null)
          })), f.$formatters.push((function(t) {
            var n, r = {};
            return $ || f.$setValidity("editable", !0), T ? (r.$model = t, T(e, r)) : (r[I.itemName] = t, n = I.viewMapper(e, r), r[I
              .itemName] = void 0, n !== I.viewMapper(e, r) ? n : t)
          }))
        }
      }
    ]).directive("uibTypeahead", (function() {
      return {
        controller: "UibTypeaheadController",
        require: ["ngModel", "uibTypeahead"],
        link: function(e, t, n, r) {
          r[1].init(r[0])
        }
      }
    })).directive("uibTypeaheadPopup", ["$$debounce", function(e) {
      return {
        scope: {
          matches: "=",
          query: "=",
          active: "=",
          position: "&",
          moveInProgress: "=",
          select: "&",
          assignIsOpen: "&",
          debounce: "&"
        },
        replace: !0,
        templateUrl: function(e, t) {
          return t.popupTemplateUrl || "uib/template/typeahead/typeahead-popup.html"
        },
        link: function(t, n, r) {
          t.templateUrl = r.templateUrl, t.isOpen = function() {
            var e = t.matches.length > 0;
            return t.assignIsOpen({
              isOpen: e
            }), e
          }, t.isActive = function(e) {
            return t.active === e
          }, t.selectActive = function(e) {
            t.active = e
          }, t.selectMatch = function(n, r) {
            var i = t.debounce();
            angular.isNumber(i) || angular.isObject(i) ? e((function() {
              t.select({
                activeIdx: n,
                evt: r
              })
            }), angular.isNumber(i) ? i : i.default) : t.select({
              activeIdx: n,
              evt: r
            })
          }
        }
      }
    }]).directive("uibTypeaheadMatch", ["$templateRequest", "$compile", "$parse", function(e, t, n) {
      return {
        scope: {
          index: "=",
          match: "=",
          query: "="
        },
        link: function(r, i, o) {
          var a = n(o.templateUrl)(r.$parent) || "uib/template/typeahead/typeahead-match.html";
          e(a).then((function(e) {
            var n = angular.element(e.trim());
            i.replaceWith(n), t(n)(r)
          }))
        }
      }
    }]).filter("uibTypeaheadHighlight", ["$sce", "$injector", "$log", function(e, t, n) {
      var r;
      return r = t.has("$sanitize"),
        function(t, i) {
          return !r && function(e) {
            return /<.*>/g.test(e)
          }(t) && n.warn("Unsafe use of typeahead please use ngSanitize"), t = i ? ("" + t).replace(new RegExp(i.replace(/([.?*+^$[\]\\(){}|-])/g,
            "\\$1"), "gi"), "<strong>$&</strong>") : t, r || (t = e.trustAsHtml(t)), t
        }
    }]), angular.module("uib/template/accordion/accordion-group.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/accordion/accordion-group.html",
        '<div role="tab" id="{{::headingId}}" aria-selected="{{isOpen}}" class="panel-heading" ng-keypress="toggleOpen($event)">\n  <h4 class="panel-title">\n    <a role="button" data-toggle="collapse" href aria-expanded="{{isOpen}}" aria-controls="{{::panelId}}" tabindex="0" class="accordion-toggle" ng-click="toggleOpen()" uib-accordion-transclude="heading" ng-disabled="isDisabled" uib-tabindex-toggle><span uib-accordion-header ng-class="{\'text-muted\': isDisabled}">{{heading}}</span></a>\n  </h4>\n</div>\n<div id="{{::panelId}}" aria-labelledby="{{::headingId}}" aria-hidden="{{!isOpen}}" role="tabpanel" class="panel-collapse collapse" uib-collapse="!isOpen">\n  <div class="panel-body" ng-transclude></div>\n</div>\n'
        )
    }]), angular.module("uib/template/accordion/accordion.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/accordion/accordion.html", '<div role="tablist" class="panel-group" ng-transclude></div>')
    }]), angular.module("uib/template/alert/alert.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/alert/alert.html",
        '<button ng-show="closeable" type="button" class="close" ng-click="close({$event: $event})">\n  <span aria-hidden="true">&times;</span>\n  <span class="sr-only">Close</span>\n</button>\n<div ng-transclude></div>\n'
        )
    }]), angular.module("uib/template/carousel/carousel.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/carousel/carousel.html",
        '<div class="carousel-inner" ng-transclude></div>\n<a role="button" href class="left carousel-control" ng-click="prev()" ng-class="{ disabled: isPrevDisabled() }" ng-show="slides.length > 1">\n  <span aria-hidden="true" class="glyphicon glyphicon-chevron-left"></span>\n  <span class="sr-only">previous</span>\n</a>\n<a role="button" href class="right carousel-control" ng-click="next()" ng-class="{ disabled: isNextDisabled() }" ng-show="slides.length > 1">\n  <span aria-hidden="true" class="glyphicon glyphicon-chevron-right"></span>\n  <span class="sr-only">next</span>\n</a>\n<ol class="carousel-indicators" ng-show="slides.length > 1">\n  <li ng-repeat="slide in slides | orderBy:indexOfSlide track by $index" ng-class="{ active: isActive(slide) }" ng-click="select(slide)">\n    <span class="sr-only">slide {{ $index + 1 }} of {{ slides.length }}<span ng-if="isActive(slide)">, currently active</span></span>\n  </li>\n</ol>\n'
        )
    }]), angular.module("uib/template/carousel/slide.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/carousel/slide.html", '<div class="text-center" ng-transclude></div>\n')
    }]), angular.module("uib/template/datepicker/datepicker.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepicker/datepicker.html",
        '<div ng-switch="datepickerMode">\n  <div uib-daypicker ng-switch-when="day" tabindex="0" class="uib-daypicker"></div>\n  <div uib-monthpicker ng-switch-when="month" tabindex="0" class="uib-monthpicker"></div>\n  <div uib-yearpicker ng-switch-when="year" tabindex="0" class="uib-yearpicker"></div>\n</div>\n'
        )
    }]), angular.module("uib/template/datepicker/day.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepicker/day.html",
        '<table role="grid" aria-labelledby="{{::uniqueId}}-title" aria-activedescendant="{{activeDateId}}">\n  <thead>\n    <tr>\n      <th><button type="button" class="btn btn-default btn-sm pull-left uib-left" ng-click="move(-1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-left"></i><span class="sr-only">previous</span></button></th>\n      <th colspan="{{::5 + showWeeks}}"><button id="{{::uniqueId}}-title" role="heading" aria-live="assertive" aria-atomic="true" type="button" class="btn btn-default btn-sm uib-title" ng-click="toggleMode()" ng-disabled="datepickerMode === maxMode" tabindex="-1"><strong>{{title}}</strong></button></th>\n      <th><button type="button" class="btn btn-default btn-sm pull-right uib-right" ng-click="move(1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-right"></i><span class="sr-only">next</span></button></th>\n    </tr>\n    <tr>\n      <th ng-if="showWeeks" class="text-center"></th>\n      <th ng-repeat="label in ::labels track by $index" class="text-center"><small aria-label="{{::label.full}}">{{::label.abbr}}</small></th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr class="uib-weeks" ng-repeat="row in rows track by $index" role="row">\n      <td ng-if="showWeeks" class="text-center h6"><em>{{ weekNumbers[$index] }}</em></td>\n      <td ng-repeat="dt in row" class="uib-day text-center" role="gridcell"\n        id="{{::dt.uid}}"\n        ng-class="::dt.customClass">\n        <button type="button" class="btn btn-default btn-sm"\n          uib-is-class="\n            \'btn-info\' for selectedDt,\n            \'active\' for activeDt\n            on dt"\n          ng-click="select(dt.date)"\n          ng-disabled="::dt.disabled"\n          tabindex="-1"><span ng-class="::{\'text-muted\': dt.secondary, \'text-info\': dt.current}">{{::dt.label}}</span></button>\n      </td>\n    </tr>\n  </tbody>\n</table>\n'
        )
    }]), angular.module("uib/template/datepicker/month.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepicker/month.html",
        '<table role="grid" aria-labelledby="{{::uniqueId}}-title" aria-activedescendant="{{activeDateId}}">\n  <thead>\n    <tr>\n      <th><button type="button" class="btn btn-default btn-sm pull-left uib-left" ng-click="move(-1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-left"></i><span class="sr-only">previous</span></button></th>\n      <th colspan="{{::yearHeaderColspan}}"><button id="{{::uniqueId}}-title" role="heading" aria-live="assertive" aria-atomic="true" type="button" class="btn btn-default btn-sm uib-title" ng-click="toggleMode()" ng-disabled="datepickerMode === maxMode" tabindex="-1"><strong>{{title}}</strong></button></th>\n      <th><button type="button" class="btn btn-default btn-sm pull-right uib-right" ng-click="move(1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-right"></i><span class="sr-only">next</span></i></button></th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr class="uib-months" ng-repeat="row in rows track by $index" role="row">\n      <td ng-repeat="dt in row" class="uib-month text-center" role="gridcell"\n        id="{{::dt.uid}}"\n        ng-class="::dt.customClass">\n        <button type="button" class="btn btn-default"\n          uib-is-class="\n            \'btn-info\' for selectedDt,\n            \'active\' for activeDt\n            on dt"\n          ng-click="select(dt.date)"\n          ng-disabled="::dt.disabled"\n          tabindex="-1"><span ng-class="::{\'text-info\': dt.current}">{{::dt.label}}</span></button>\n      </td>\n    </tr>\n  </tbody>\n</table>\n'
        )
    }]), angular.module("uib/template/datepicker/popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepicker/popup.html",
        '<div>\n  <ul class="uib-datepicker-popup dropdown-menu uib-position-measure" dropdown-nested ng-if="isOpen" ng-keydown="keydown($event)" ng-click="$event.stopPropagation()">\n    <li ng-transclude></li>\n    <li ng-if="showButtonBar" class="uib-button-bar">\n      <span class="btn-group pull-left">\n        <button type="button" class="btn btn-sm btn-info uib-datepicker-current" ng-click="select(\'today\', $event)" ng-disabled="isDisabled(\'today\')">{{ getText(\'current\') }}</button>\n        <button type="button" class="btn btn-sm btn-danger uib-clear" ng-click="select(null, $event)">{{ getText(\'clear\') }}</button>\n      </span>\n      <button type="button" class="btn btn-sm btn-success pull-right uib-close" ng-click="close($event)">{{ getText(\'close\') }}</button>\n    </li>\n  </ul>\n</div>\n'
        )
    }]), angular.module("uib/template/datepicker/year.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepicker/year.html",
        '<table role="grid" aria-labelledby="{{::uniqueId}}-title" aria-activedescendant="{{activeDateId}}">\n  <thead>\n    <tr>\n      <th><button type="button" class="btn btn-default btn-sm pull-left uib-left" ng-click="move(-1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-left"></i><span class="sr-only">previous</span></button></th>\n      <th colspan="{{::columns - 2}}"><button id="{{::uniqueId}}-title" role="heading" aria-live="assertive" aria-atomic="true" type="button" class="btn btn-default btn-sm uib-title" ng-click="toggleMode()" ng-disabled="datepickerMode === maxMode" tabindex="-1"><strong>{{title}}</strong></button></th>\n      <th><button type="button" class="btn btn-default btn-sm pull-right uib-right" ng-click="move(1)" tabindex="-1"><i aria-hidden="true" class="glyphicon glyphicon-chevron-right"></i><span class="sr-only">next</span></button></th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr class="uib-years" ng-repeat="row in rows track by $index" role="row">\n      <td ng-repeat="dt in row" class="uib-year text-center" role="gridcell"\n        id="{{::dt.uid}}"\n        ng-class="::dt.customClass">\n        <button type="button" class="btn btn-default"\n          uib-is-class="\n            \'btn-info\' for selectedDt,\n            \'active\' for activeDt\n            on dt"\n          ng-click="select(dt.date)"\n          ng-disabled="::dt.disabled"\n          tabindex="-1"><span ng-class="::{\'text-info\': dt.current}">{{::dt.label}}</span></button>\n      </td>\n    </tr>\n  </tbody>\n</table>\n'
        )
    }]), angular.module("uib/template/datepickerPopup/popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/datepickerPopup/popup.html",
        '<ul role="presentation" class="uib-datepicker-popup dropdown-menu uib-position-measure" dropdown-nested ng-if="isOpen" ng-keydown="keydown($event)" ng-click="$event.stopPropagation()">\n  <li ng-transclude></li>\n  <li ng-if="showButtonBar" class="uib-button-bar">\n    <span class="btn-group pull-left">\n      <button type="button" class="btn btn-sm btn-info uib-datepicker-current" ng-click="select(\'today\', $event)" ng-disabled="isDisabled(\'today\')">{{ getText(\'current\') }}</button>\n      <button type="button" class="btn btn-sm btn-danger uib-clear" ng-click="select(null, $event)">{{ getText(\'clear\') }}</button>\n    </span>\n    <button type="button" class="btn btn-sm btn-success pull-right uib-close" ng-click="close($event)">{{ getText(\'close\') }}</button>\n  </li>\n</ul>\n'
        )
    }]), angular.module("uib/template/modal/backdrop.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/modal/backdrop.html",
        '<div class="modal-backdrop"\n     uib-modal-animation-class="fade"\n     modal-in-class="in"\n     ng-style="{\'z-index\': 1040 + (index && 1 || 0) + index*10}"\n></div>\n'
        )
    }]), angular.module("uib/template/modal/window.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/modal/window.html",
        "<div class=\"modal-dialog {{size ? 'modal-' + size : ''}}\"><div class=\"modal-content\" uib-modal-transclude></div></div>\n")
    }]), angular.module("uib/template/pager/pager.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/pager/pager.html",
        '<li ng-class="{disabled: noPrevious()||ngDisabled, previous: align}"><a href ng-click="selectPage(page - 1, $event)" ng-disabled="noPrevious()||ngDisabled" uib-tabindex-toggle>{{::getText(\'previous\')}}</a></li>\n<li ng-class="{disabled: noNext()||ngDisabled, next: align}"><a href ng-click="selectPage(page + 1, $event)" ng-disabled="noNext()||ngDisabled" uib-tabindex-toggle>{{::getText(\'next\')}}</a></li>\n'
        )
    }]), angular.module("uib/template/pagination/pagination.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/pagination/pagination.html",
        '<li role="menuitem" ng-if="::boundaryLinks" ng-class="{disabled: noPrevious()||ngDisabled}" class="pagination-first"><a href ng-click="selectPage(1, $event)" ng-disabled="noPrevious()||ngDisabled" uib-tabindex-toggle>{{::getText(\'first\')}}</a></li>\n<li role="menuitem" ng-if="::directionLinks" ng-class="{disabled: noPrevious()||ngDisabled}" class="pagination-prev"><a href ng-click="selectPage(page - 1, $event)" ng-disabled="noPrevious()||ngDisabled" uib-tabindex-toggle>{{::getText(\'previous\')}}</a></li>\n<li role="menuitem" ng-repeat="page in pages track by $index" ng-class="{active: page.active,disabled: ngDisabled&&!page.active}" class="pagination-page"><a href ng-click="selectPage(page.number, $event)" ng-disabled="ngDisabled&&!page.active" uib-tabindex-toggle>{{page.text}}</a></li>\n<li role="menuitem" ng-if="::directionLinks" ng-class="{disabled: noNext()||ngDisabled}" class="pagination-next"><a href ng-click="selectPage(page + 1, $event)" ng-disabled="noNext()||ngDisabled" uib-tabindex-toggle>{{::getText(\'next\')}}</a></li>\n<li role="menuitem" ng-if="::boundaryLinks" ng-class="{disabled: noNext()||ngDisabled}" class="pagination-last"><a href ng-click="selectPage(totalPages, $event)" ng-disabled="noNext()||ngDisabled" uib-tabindex-toggle>{{::getText(\'last\')}}</a></li>\n'
        )
    }]), angular.module("uib/template/tooltip/tooltip-html-popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/tooltip/tooltip-html-popup.html",
        '<div class="tooltip-arrow"></div>\n<div class="tooltip-inner" ng-bind-html="contentExp()"></div>\n')
    }]), angular.module("uib/template/tooltip/tooltip-popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/tooltip/tooltip-popup.html", '<div class="tooltip-arrow"></div>\n<div class="tooltip-inner" ng-bind="content"></div>\n')
    }]), angular.module("uib/template/tooltip/tooltip-template-popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/tooltip/tooltip-template-popup.html",
        '<div class="tooltip-arrow"></div>\n<div class="tooltip-inner"\n  uib-tooltip-template-transclude="contentExp()"\n  tooltip-template-transclude-scope="originScope()"></div>\n'
        )
    }]), angular.module("uib/template/popover/popover-html.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/popover/popover-html.html",
        '<div class="arrow"></div>\n\n<div class="popover-inner">\n    <h3 class="popover-title" ng-bind="uibTitle" ng-if="uibTitle"></h3>\n    <div class="popover-content" ng-bind-html="contentExp()"></div>\n</div>\n'
        )
    }]), angular.module("uib/template/popover/popover-template.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/popover/popover-template.html",
        '<div class="arrow"></div>\n\n<div class="popover-inner">\n    <h3 class="popover-title" ng-bind="uibTitle" ng-if="uibTitle"></h3>\n    <div class="popover-content"\n      uib-tooltip-template-transclude="contentExp()"\n      tooltip-template-transclude-scope="originScope()"></div>\n</div>\n'
        )
    }]), angular.module("uib/template/popover/popover.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/popover/popover.html",
        '<div class="arrow"></div>\n\n<div class="popover-inner">\n    <h3 class="popover-title" ng-bind="uibTitle" ng-if="uibTitle"></h3>\n    <div class="popover-content" ng-bind="content"></div>\n</div>\n'
        )
    }]), angular.module("uib/template/progressbar/bar.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/progressbar/bar.html",
        '<div class="progress-bar" ng-class="type && \'progress-bar-\' + type" role="progressbar" aria-valuenow="{{value}}" aria-valuemin="0" aria-valuemax="{{max}}" ng-style="{width: (percent < 100 ? percent : 100) + \'%\'}" aria-valuetext="{{percent | number:0}}%" aria-labelledby="{{::title}}" ng-transclude></div>\n'
        )
    }]), angular.module("uib/template/progressbar/progress.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/progressbar/progress.html", '<div class="progress" ng-transclude aria-labelledby="{{::title}}"></div>')
    }]), angular.module("uib/template/progressbar/progressbar.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/progressbar/progressbar.html",
        '<div class="progress">\n  <div class="progress-bar" ng-class="type && \'progress-bar-\' + type" role="progressbar" aria-valuenow="{{value}}" aria-valuemin="0" aria-valuemax="{{max}}" ng-style="{width: (percent < 100 ? percent : 100) + \'%\'}" aria-valuetext="{{percent | number:0}}%" aria-labelledby="{{::title}}" ng-transclude></div>\n</div>\n'
        )
    }]), angular.module("uib/template/rating/rating.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/rating/rating.html",
        '<span ng-mouseleave="reset()" ng-keydown="onKeydown($event)" tabindex="0" role="slider" aria-valuemin="0" aria-valuemax="{{range.length}}" aria-valuenow="{{value}}" aria-valuetext="{{title}}">\n    <span ng-repeat-start="r in range track by $index" class="sr-only">({{ $index < value ? \'*\' : \' \' }})</span>\n    <i ng-repeat-end ng-mouseenter="enter($index + 1)" ng-click="rate($index + 1)" class="glyphicon" ng-class="$index < value && (r.stateOn || \'glyphicon-star\') || (r.stateOff || \'glyphicon-star-empty\')" ng-attr-title="{{r.title}}"></i>\n</span>\n'
        )
    }]), angular.module("uib/template/tabs/tab.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/tabs/tab.html",
        '<li ng-class="[{active: active, disabled: disabled}, classes]" class="uib-tab nav-item">\n  <a href ng-click="select($event)" class="nav-link" uib-tab-heading-transclude>{{heading}}</a>\n</li>\n'
        )
    }]), angular.module("uib/template/tabs/tabset.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/tabs/tabset.html",
        '<div>\n  <ul class="nav nav-{{tabset.type || \'tabs\'}}" ng-class="{\'nav-stacked\': vertical, \'nav-justified\': justified}" ng-transclude></ul>\n  <div class="tab-content">\n    <div class="tab-pane"\n         ng-repeat="tab in tabset.tabs"\n         ng-class="{active: tabset.active === tab.index}"\n         uib-tab-content-transclude="tab">\n    </div>\n  </div>\n</div>\n'
        )
    }]), angular.module("uib/template/timepicker/timepicker.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/timepicker/timepicker.html",
        '<table class="uib-timepicker">\n  <tbody>\n    <tr class="text-center" ng-show="::showSpinners">\n      <td class="uib-increment hours"><a ng-click="incrementHours()" ng-class="{disabled: noIncrementHours()}" class="btn btn-link" ng-disabled="noIncrementHours()" tabindex="-1"><span class="glyphicon glyphicon-chevron-up"></span></a></td>\n      <td>&nbsp;</td>\n      <td class="uib-increment minutes"><a ng-click="incrementMinutes()" ng-class="{disabled: noIncrementMinutes()}" class="btn btn-link" ng-disabled="noIncrementMinutes()" tabindex="-1"><span class="glyphicon glyphicon-chevron-up"></span></a></td>\n      <td ng-show="showSeconds">&nbsp;</td>\n      <td ng-show="showSeconds" class="uib-increment seconds"><a ng-click="incrementSeconds()" ng-class="{disabled: noIncrementSeconds()}" class="btn btn-link" ng-disabled="noIncrementSeconds()" tabindex="-1"><span class="glyphicon glyphicon-chevron-up"></span></a></td>\n      <td ng-show="showMeridian"></td>\n    </tr>\n    <tr>\n      <td class="form-group uib-time hours" ng-class="{\'has-error\': invalidHours}">\n        <input type="text" placeholder="HH" ng-model="hours" ng-change="updateHours()" class="form-control text-center" ng-readonly="::readonlyInput" maxlength="2" tabindex="{{::tabindex}}" ng-disabled="noIncrementHours()" ng-blur="blur()">\n      </td>\n      <td class="uib-separator">:</td>\n      <td class="form-group uib-time minutes" ng-class="{\'has-error\': invalidMinutes}">\n        <input type="text" placeholder="MM" ng-model="minutes" ng-change="updateMinutes()" class="form-control text-center" ng-readonly="::readonlyInput" maxlength="2" tabindex="{{::tabindex}}" ng-disabled="noIncrementMinutes()" ng-blur="blur()">\n      </td>\n      <td ng-show="showSeconds" class="uib-separator">:</td>\n      <td class="form-group uib-time seconds" ng-class="{\'has-error\': invalidSeconds}" ng-show="showSeconds">\n        <input type="text" placeholder="SS" ng-model="seconds" ng-change="updateSeconds()" class="form-control text-center" ng-readonly="readonlyInput" maxlength="2" tabindex="{{::tabindex}}" ng-disabled="noIncrementSeconds()" ng-blur="blur()">\n      </td>\n      <td ng-show="showMeridian" class="uib-time am-pm"><button type="button" ng-class="{disabled: noToggleMeridian()}" class="btn btn-default text-center" ng-click="toggleMeridian()" ng-disabled="noToggleMeridian()" tabindex="{{::tabindex}}">{{meridian}}</button></td>\n    </tr>\n    <tr class="text-center" ng-show="::showSpinners">\n      <td class="uib-decrement hours"><a ng-click="decrementHours()" ng-class="{disabled: noDecrementHours()}" class="btn btn-link" ng-disabled="noDecrementHours()" tabindex="-1"><span class="glyphicon glyphicon-chevron-down"></span></a></td>\n      <td>&nbsp;</td>\n      <td class="uib-decrement minutes"><a ng-click="decrementMinutes()" ng-class="{disabled: noDecrementMinutes()}" class="btn btn-link" ng-disabled="noDecrementMinutes()" tabindex="-1"><span class="glyphicon glyphicon-chevron-down"></span></a></td>\n      <td ng-show="showSeconds">&nbsp;</td>\n      <td ng-show="showSeconds" class="uib-decrement seconds"><a ng-click="decrementSeconds()" ng-class="{disabled: noDecrementSeconds()}" class="btn btn-link" ng-disabled="noDecrementSeconds()" tabindex="-1"><span class="glyphicon glyphicon-chevron-down"></span></a></td>\n      <td ng-show="showMeridian"></td>\n    </tr>\n  </tbody>\n</table>\n'
        )
    }]), angular.module("uib/template/typeahead/typeahead-match.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/typeahead/typeahead-match.html",
        '<a href\n   tabindex="-1"\n   ng-bind-html="match.label | uibTypeaheadHighlight:query"\n   ng-attr-title="{{match.label}}"></a>\n')
    }]), angular.module("uib/template/typeahead/typeahead-popup.html", []).run(["$templateCache", function(e) {
      e.put("uib/template/typeahead/typeahead-popup.html",
        '<ul class="dropdown-menu" ng-show="isOpen() && !moveInProgress" ng-style="{top: position().top+\'px\', left: position().left+\'px\'}" role="listbox" aria-hidden="{{!isOpen()}}">\n    <li class="uib-typeahead-match" ng-repeat="match in matches track by $index" ng-class="{active: isActive($index) }" ng-mouseenter="selectActive($index)" ng-click="selectMatch($index, $event)" role="option" id="{{::match.id}}">\n        <div uib-typeahead-match index="$index" match="match" query="query" template-url="templateUrl"></div>\n    </li>\n</ul>\n'
        )
    }]), angular.module("ui.bootstrap.carousel").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibCarouselCss && angular.element(document).find("head").prepend(
        '<style type="text/css">.ng-animate.item:not(.left):not(.right){-webkit-transition:0s ease-in-out left;transition:0s ease-in-out left}</style>'
        ), angular.$$uibCarouselCss = !0
    })), angular.module("ui.bootstrap.datepicker").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibDatepickerCss && angular.element(document).find("head").prepend(
        '<style type="text/css">.uib-datepicker .uib-title{width:100%;}.uib-day button,.uib-month button,.uib-year button{min-width:100%;}.uib-left,.uib-right{width:100%}</style>'
        ), angular.$$uibDatepickerCss = !0
    })), angular.module("ui.bootstrap.position").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibPositionCss && angular.element(document).find("head").prepend(
        '<style type="text/css">.uib-position-measure{display:block !important;visibility:hidden !important;position:absolute !important;top:-9999px !important;left:-9999px !important;}.uib-position-scrollbar-measure{position:absolute !important;top:-9999px !important;width:50px !important;height:50px !important;overflow:scroll !important;}.uib-position-body-scrollbar-measure{overflow:scroll !important;}</style>'
        ), angular.$$uibPositionCss = !0
    })), angular.module("ui.bootstrap.datepickerPopup").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibDatepickerpopupCss && angular.element(document).find("head").prepend(
        '<style type="text/css">.uib-datepicker-popup.dropdown-menu{display:block;float:none;margin:0;}.uib-button-bar{padding:10px 9px 2px;}</style>'
        ), angular.$$uibDatepickerpopupCss = !0
    })), angular.module("ui.bootstrap.tooltip").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibTooltipCss && angular.element(document).find("head").prepend(
        '<style type="text/css">[uib-tooltip-popup].tooltip.top-left > .tooltip-arrow,[uib-tooltip-popup].tooltip.top-right > .tooltip-arrow,[uib-tooltip-popup].tooltip.bottom-left > .tooltip-arrow,[uib-tooltip-popup].tooltip.bottom-right > .tooltip-arrow,[uib-tooltip-popup].tooltip.left-top > .tooltip-arrow,[uib-tooltip-popup].tooltip.left-bottom > .tooltip-arrow,[uib-tooltip-popup].tooltip.right-top > .tooltip-arrow,[uib-tooltip-popup].tooltip.right-bottom > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.top-left > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.top-right > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.bottom-left > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.bottom-right > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.left-top > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.left-bottom > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.right-top > .tooltip-arrow,[uib-tooltip-html-popup].tooltip.right-bottom > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.top-left > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.top-right > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.bottom-left > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.bottom-right > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.left-top > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.left-bottom > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.right-top > .tooltip-arrow,[uib-tooltip-template-popup].tooltip.right-bottom > .tooltip-arrow,[uib-popover-popup].popover.top-left > .arrow,[uib-popover-popup].popover.top-right > .arrow,[uib-popover-popup].popover.bottom-left > .arrow,[uib-popover-popup].popover.bottom-right > .arrow,[uib-popover-popup].popover.left-top > .arrow,[uib-popover-popup].popover.left-bottom > .arrow,[uib-popover-popup].popover.right-top > .arrow,[uib-popover-popup].popover.right-bottom > .arrow,[uib-popover-html-popup].popover.top-left > .arrow,[uib-popover-html-popup].popover.top-right > .arrow,[uib-popover-html-popup].popover.bottom-left > .arrow,[uib-popover-html-popup].popover.bottom-right > .arrow,[uib-popover-html-popup].popover.left-top > .arrow,[uib-popover-html-popup].popover.left-bottom > .arrow,[uib-popover-html-popup].popover.right-top > .arrow,[uib-popover-html-popup].popover.right-bottom > .arrow,[uib-popover-template-popup].popover.top-left > .arrow,[uib-popover-template-popup].popover.top-right > .arrow,[uib-popover-template-popup].popover.bottom-left > .arrow,[uib-popover-template-popup].popover.bottom-right > .arrow,[uib-popover-template-popup].popover.left-top > .arrow,[uib-popover-template-popup].popover.left-bottom > .arrow,[uib-popover-template-popup].popover.right-top > .arrow,[uib-popover-template-popup].popover.right-bottom > .arrow{top:auto;bottom:auto;left:auto;right:auto;margin:0;}[uib-popover-popup].popover,[uib-popover-html-popup].popover,[uib-popover-template-popup].popover{display:block !important;}</style>'
        ), angular.$$uibTooltipCss = !0
    })), angular.module("ui.bootstrap.timepicker").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibTimepickerCss && angular.element(document).find("head").prepend(
        '<style type="text/css">.uib-time input{width:50px;}</style>'), angular.$$uibTimepickerCss = !0
    })), angular.module("ui.bootstrap.typeahead").run((function() {
      !angular.$$csp().noInlineStyle && !angular.$$uibTypeaheadCss && angular.element(document).find("head").prepend(
        '<style type="text/css">[uib-typeahead-popup].dropdown-menu{display:block;}</style>'), angular.$$uibTypeaheadCss = !0
    }))
  }, {}],
  30: [function(e, t, n) {
    e("./dist/ui-bootstrap-tpls"), t.exports = "ui.bootstrap"
  }, {
    "./dist/ui-bootstrap-tpls": 29
  }],
  31: [function(e, t, n) {
    /**
     * @license AngularJS v1.8.2
     * (c) 2010-2020 Google LLC. http://angularjs.org
     * License: MIT
     */
    ! function(e) {
      "use strict";
      var t = {
        objectMaxDepth: 5,
        urlErrorParamsEnabled: !0
      };

      function n(e) {
        if (!V(e)) return t;
        _(e.objectMaxDepth) && (t.objectMaxDepth = r(e.objectMaxDepth) ? e.objectMaxDepth : NaN), _(e.urlErrorParamsEnabled) && X(e
          .urlErrorParamsEnabled) && (t.urlErrorParamsEnabled = e.urlErrorParamsEnabled)
      }

      function r(e) {
        return F(e) && e > 0
      }

      function i(e, n) {
        n = n || Error;
        var r = "https://errors.angularjs.org/1.8.2/",
          i = r.replace(".", "\\.") + "[\\s\\S]*",
          o = new RegExp(i, "g");
        return function() {
          var i, a, s = arguments[0],
            u = arguments[1],
            l = "[" + (e ? e + ":" : "") + s + "] ",
            c = pe(arguments, 2).map((function(e) {
              return Be(e, t.objectMaxDepth)
            }));
          if (l += u.replace(/\{\d+\}/g, (function(e) {
              var t = +e.slice(1, -1);
              return t < c.length ? c[t].replace(o, "") : e
            })), l += "\n" + r + (e ? e + "/" : "") + s, t.urlErrorParamsEnabled)
            for (a = 0, i = "?"; a < c.length; a++, i = "&") l += i + "p" + a + "=" + encodeURIComponent(c[a]);
          return new n(l)
        }
      }
      var o, a, s, u, l = /^\/(.+)\/([a-z]*)$/,
        c = "validity",
        d = Object.prototype.hasOwnProperty,
        p = function(e) {
          return H(e) ? e.toLowerCase() : e
        },
        f = function(e) {
          return H(e) ? e.toUpperCase() : e
        },
        h = [].slice,
        g = [].splice,
        m = [].push,
        v = Object.prototype.toString,
        $ = Object.getPrototypeOf,
        b = i("ng"),
        y = e.angular || (e.angular = {}),
        w = 0;

      function x(e) {
        if (null == e || K(e)) return !1;
        if (G(e) || H(e) || a && e instanceof a) return !0;
        var t = "length" in Object(e) && e.length;
        return F(t) && (t >= 0 && t - 1 in e || "function" == typeof e.item)
      }

      function C(e, t, n) {
        var r, i;
        if (e)
          if (z(e))
            for (r in e) "prototype" !== r && "length" !== r && "name" !== r && e.hasOwnProperty(r) && t.call(n, e[r], r, e);
          else if (G(e) || x(e)) {
          var o = "object" != typeof e;
          for (r = 0, i = e.length; r < i; r++)(o || r in e) && t.call(n, e[r], r, e)
        } else if (e.forEach && e.forEach !== C) e.forEach(t, n, e);
        else if (U(e))
          for (r in e) t.call(n, e[r], r, e);
        else if ("function" == typeof e.hasOwnProperty)
          for (r in e) e.hasOwnProperty(r) && t.call(n, e[r], r, e);
        else
          for (r in e) d.call(e, r) && t.call(n, e[r], r, e);
        return e
      }

      function S(e, t, n) {
        for (var r = Object.keys(e).sort(), i = 0; i < r.length; i++) t.call(n, e[r[i]], r[i]);
        return r
      }

      function A(e) {
        return function(t, n) {
          e(n, t)
        }
      }

      function T() {
        return ++w
      }

      function k(e, t) {
        t ? e.$$hashKey = t : delete e.$$hashKey
      }

      function D(e, t, n) {
        for (var r = e.$$hashKey, i = 0, o = t.length; i < o; ++i) {
          var a = t[i];
          if (V(a) || z(a))
            for (var s = Object.keys(a), u = 0, l = s.length; u < l; u++) {
              var c = s[u],
                d = a[c];
              n && V(d) ? B(d) ? e[c] = new Date(d.valueOf()) : Y(d) ? e[c] = new RegExp(d) : d.nodeName ? e[c] = d.cloneNode(!0) : ne(d) ? e[c] = d
                .clone() : "__proto__" !== c && (V(e[c]) || (e[c] = G(d) ? [] : {}), D(e[c], [d], !0)) : e[c] = d
            }
        }
        return k(e, r), e
      }

      function M(e) {
        return D(e, h.call(arguments, 1), !1)
      }

      function E(e) {
        return D(e, h.call(arguments, 1), !0)
      }

      function O(e) {
        return parseInt(e, 10)
      }
      o = e.document.documentMode;
      var N = Number.isNaN || function(e) {
        return e != e
      };

      function P(e, t) {
        return M(Object.create(e), t)
      }

      function q() {}

      function I(e) {
        return e
      }

      function L(e) {
        return function() {
          return e
        }
      }

      function R(e) {
        return z(e.toString) && e.toString !== v
      }

      function j(e) {
        return void 0 === e
      }

      function _(e) {
        return void 0 !== e
      }

      function V(e) {
        return null !== e && "object" == typeof e
      }

      function U(e) {
        return null !== e && "object" == typeof e && !$(e)
      }

      function H(e) {
        return "string" == typeof e
      }

      function F(e) {
        return "number" == typeof e
      }

      function B(e) {
        return "[object Date]" === v.call(e)
      }

      function G(e) {
        return Array.isArray(e) || e instanceof Array
      }

      function W(e) {
        switch (v.call(e)) {
          case "[object Error]":
          case "[object Exception]":
          case "[object DOMException]":
            return !0;
          default:
            return e instanceof Error
        }
      }

      function z(e) {
        return "function" == typeof e
      }

      function Y(e) {
        return "[object RegExp]" === v.call(e)
      }

      function K(e) {
        return e && e.window === e
      }

      function Q(e) {
        return e && e.$evalAsync && e.$watch
      }

      function X(e) {
        return "boolean" == typeof e
      }

      function J(e) {
        return e && z(e.then)
      }
      q.$inject = [], I.$inject = [];
      var Z = /^\[object (?:Uint8|Uint8Clamped|Uint16|Uint32|Int8|Int16|Int32|Float32|Float64)Array]$/;
      var ee = function(e) {
          return H(e) ? e.trim() : e
        },
        te = function(e) {
          return e.replace(/([-()[\]{}+?*.$^|,:#<!\\])/g, "\\$1").replace(/\x08/g, "\\x08")
        };

      function ne(e) {
        return !(!e || !(e.nodeName || e.prop && e.attr && e.find))
      }

      function re(e) {
        return p(e.nodeName || e[0] && e[0].nodeName)
      }

      function ie(e, t) {
        return -1 !== Array.prototype.indexOf.call(e, t)
      }

      function oe(e, t) {
        var n = e.indexOf(t);
        return n >= 0 && e.splice(n, 1), n
      }

      function ae(e, t, n) {
        var i, o, a = [],
          s = [];
        if (n = r(n) ? n : NaN, t) {
          if ((o = t) && F(o.length) && Z.test(v.call(o)) || (i = t, "[object ArrayBuffer]" === v.call(i))) throw b("cpta",
            "Can't copy! TypedArray destination cannot be mutated.");
          if (e === t) throw b("cpi", "Can't copy! Source and destination are identical.");
          return G(t) ? t.length = 0 : C(t, (function(e, n) {
            "$$hashKey" !== n && delete t[n]
          })), a.push(e), s.push(t), u(e, t, n)
        }
        return l(e, n);

        function u(e, t, n) {
          if (--n < 0) return "...";
          var r, i = t.$$hashKey;
          if (G(e))
            for (var o = 0, a = e.length; o < a; o++) t.push(l(e[o], n));
          else if (U(e))
            for (r in e) t[r] = l(e[r], n);
          else if (e && "function" == typeof e.hasOwnProperty)
            for (r in e) e.hasOwnProperty(r) && (t[r] = l(e[r], n));
          else
            for (r in e) d.call(e, r) && (t[r] = l(e[r], n));
          return k(t, i), t
        }

        function l(e, t) {
          if (!V(e)) return e;
          var n = a.indexOf(e);
          if (-1 !== n) return s[n];
          if (K(e) || Q(e)) throw b("cpws", "Can't copy! Making copies of Window or Scope instances is not supported.");
          var r = !1,
            i = function(e) {
              switch (v.call(e)) {
                case "[object Int8Array]":
                case "[object Int16Array]":
                case "[object Int32Array]":
                case "[object Float32Array]":
                case "[object Float64Array]":
                case "[object Uint8Array]":
                case "[object Uint8ClampedArray]":
                case "[object Uint16Array]":
                case "[object Uint32Array]":
                  return new e.constructor(l(e.buffer), e.byteOffset, e.length);
                case "[object ArrayBuffer]":
                  if (!e.slice) {
                    var t = new ArrayBuffer(e.byteLength);
                    return new Uint8Array(t).set(new Uint8Array(e)), t
                  }
                  return e.slice(0);
                case "[object Boolean]":
                case "[object Number]":
                case "[object String]":
                case "[object Date]":
                  return new e.constructor(e.valueOf());
                case "[object RegExp]":
                  var n = new RegExp(e.source, e.toString().match(/[^/]*$/)[0]);
                  return n.lastIndex = e.lastIndex, n;
                case "[object Blob]":
                  return new e.constructor([e], {
                    type: e.type
                  })
              }
              if (z(e.cloneNode)) return e.cloneNode(!0)
            }(e);
          return void 0 === i && (i = G(e) ? [] : Object.create($(e)), r = !0), a.push(e), s.push(i), r ? u(e, i, t) : i
        }
      }

      function se(e, t) {
        return e === t || e != e && t != t
      }

      function ue(e, t) {
        if (e === t) return !0;
        if (null === e || null === t) return !1;
        if (e != e && t != t) return !0;
        var n, r, i, o = typeof e;
        if (o === typeof t && "object" === o) {
          if (!G(e)) {
            if (B(e)) return !!B(t) && se(e.getTime(), t.getTime());
            if (Y(e)) return !!Y(t) && e.toString() === t.toString();
            if (Q(e) || Q(t) || K(e) || K(t) || G(t) || B(t) || Y(t)) return !1;
            for (r in i = Ve(), e)
              if ("$" !== r.charAt(0) && !z(e[r])) {
                if (!ue(e[r], t[r])) return !1;
                i[r] = !0
              } for (r in t)
              if (!(r in i) && "$" !== r.charAt(0) && _(t[r]) && !z(t[r])) return !1;
            return !0
          }
          if (!G(t)) return !1;
          if ((n = e.length) === t.length) {
            for (r = 0; r < n; r++)
              if (!ue(e[r], t[r])) return !1;
            return !0
          }
        }
        return !1
      }
      var le = function() {
          if (!_(le.rules)) {
            var t = e.document.querySelector("[ng-csp]") || e.document.querySelector("[data-ng-csp]");
            if (t) {
              var n = t.getAttribute("ng-csp") || t.getAttribute("data-ng-csp");
              le.rules = {
                noUnsafeEval: !n || -1 !== n.indexOf("no-unsafe-eval"),
                noInlineStyle: !n || -1 !== n.indexOf("no-inline-style")
              }
            } else le.rules = {
              noUnsafeEval: function() {
                try {
                  return new Function(""), !1
                } catch (e) {
                  return !0
                }
              }(),
              noInlineStyle: !1
            }
          }
          return le.rules
        },
        ce = function() {
          if (_(ce.name_)) return ce.name_;
          var t, n, r, i, o = Te.length;
          for (n = 0; n < o; ++n)
            if (r = Te[n], t = e.document.querySelector("[" + r.replace(":", "\\:") + "jq]")) {
              i = t.getAttribute(r + "jq");
              break
            } return ce.name_ = i
        };

      function de(e, t, n) {
        return e.concat(h.call(t, n))
      }

      function pe(e, t) {
        return h.call(e, t || 0)
      }

      function fe(e, t) {
        var n = arguments.length > 2 ? pe(arguments, 2) : [];
        return !z(t) || t instanceof RegExp ? t : n.length ? function() {
          return arguments.length ? t.apply(e, de(n, arguments, 0)) : t.apply(e, n)
        } : function() {
          return arguments.length ? t.apply(e, arguments) : t.call(e)
        }
      }

      function he(t, n) {
        var r = n;
        return "string" == typeof t && "$" === t.charAt(0) && "$" === t.charAt(1) ? r = void 0 : K(n) ? r = "$WINDOW" : n && e.document === n ? r =
          "$DOCUMENT" : Q(n) && (r = "$SCOPE"), r
      }

      function ge(e, t) {
        if (!j(e)) return F(t) || (t = t ? 2 : null), JSON.stringify(e, he, t)
      }

      function me(e) {
        return H(e) ? JSON.parse(e) : e
      }
      var ve = /:/g;

      function $e(e, t) {
        e = e.replace(ve, "");
        var n = Date.parse("Jan 01, 1970 00:00:00 " + e) / 6e4;
        return N(n) ? t : n
      }

      function be(e, t) {
        return (e = new Date(e.getTime())).setMinutes(e.getMinutes() + t), e
      }

      function ye(e, t, n) {
        n = n ? -1 : 1;
        var r = e.getTimezoneOffset();
        return be(e, n * ($e(t, r) - r))
      }

      function we(e) {
        e = a(e).clone().empty();
        var t = a("<div></div>").append(e).html();
        try {
          return e[0].nodeType === He ? p(t) : t.match(/^(<[^>]+>)/)[1].replace(/^<([\w-]+)/, (function(e, t) {
            return "<" + p(t)
          }))
        } catch (e) {
          return p(t)
        }
      }

      function xe(e) {
        try {
          return decodeURIComponent(e)
        } catch (e) {}
      }

      function Ce(e) {
        var t = {};
        return C((e || "").split("&"), (function(e) {
          var n, r, i;
          e && (r = e = e.replace(/\+/g, "%20"), -1 !== (n = e.indexOf("=")) && (r = e.substring(0, n), i = e.substring(n + 1)), _(r = xe(r)) && (
            i = !_(i) || xe(i), d.call(t, r) ? G(t[r]) ? t[r].push(i) : t[r] = [t[r], i] : t[r] = i))
        })), t
      }

      function Se(e) {
        return Ae(e, !0).replace(/%26/gi, "&").replace(/%3D/gi, "=").replace(/%2B/gi, "+")
      }

      function Ae(e, t) {
        return encodeURIComponent(e).replace(/%40/gi, "@").replace(/%3A/gi, ":").replace(/%24/g, "$").replace(/%2C/gi, ",").replace(/%3B/gi, ";").replace(
          /%20/g, t ? "%20" : "+")
      }
      var Te = ["ng-", "data-ng-", "ng:", "x-ng-"];
      var ke = function(t) {
        var n = t.currentScript;
        if (!n) return !0;
        if (!(n instanceof e.HTMLScriptElement || n instanceof e.SVGScriptElement)) return !1;
        var r = n.attributes;
        return [r.getNamedItem("src"), r.getNamedItem("href"), r.getNamedItem("xlink:href")].every((function(e) {
          if (!e) return !0;
          if (!e.value) return !1;
          var n = t.createElement("a");
          if (n.href = e.value, t.location.origin === n.origin) return !0;
          switch (n.protocol) {
            case "http:":
            case "https:":
            case "ftp:":
            case "blob:":
            case "file:":
            case "data:":
              return !0;
            default:
              return !1
          }
        }))
      }(e.document);

      function De(t, n) {
        var r, i, o = {};
        if (C(Te, (function(e) {
            var n = e + "app";
            !r && t.hasAttribute && t.hasAttribute(n) && (r = t, i = t.getAttribute(n))
          })), C(Te, (function(e) {
            var n, o = e + "app";
            !r && (n = t.querySelector("[" + o.replace(":", "\\:") + "]")) && (r = n, i = n.getAttribute(o))
          })), r) {
          if (!ke) return void e.console.error(
            "AngularJS: disabling automatic bootstrap. <script> protocol indicates an extension, document.location.href does not match.");
          o.strictDi = null !== function(e, t) {
            var n, r, i = Te.length;
            for (r = 0; r < i; ++r)
              if (n = Te[r] + t, H(n = e.getAttribute(n))) return n;
            return null
          }(r, "strict-di"), n(r, i ? [i] : [], o)
        }
      }

      function Me(t, n, r) {
        V(r) || (r = {});
        r = M({
          strictDi: !1
        }, r);
        var i = function() {
            if ((t = a(t)).injector()) {
              var i = t[0] === e.document ? "document" : we(t);
              throw b("btstrpd", "App already bootstrapped with this element '{0}'", i.replace(/</, "&lt;").replace(/>/, "&gt;"))
            }(n = n || []).unshift(["$provide", function(e) {
              e.value("$rootElement", t)
            }]), r.debugInfoEnabled && n.push(["$compileProvider", function(e) {
              e.debugInfoEnabled(!0)
            }]), n.unshift("ng");
            var o = Jt(n, r.strictDi);
            return o.invoke(["$rootScope", "$rootElement", "$compile", "$injector", function(e, t, n, r) {
              e.$apply((function() {
                t.data("$injector", r), n(t)(e)
              }))
            }]), o
          },
          o = /^NG_ENABLE_DEBUG_INFO!/,
          s = /^NG_DEFER_BOOTSTRAP!/;
        if (e && o.test(e.name) && (r.debugInfoEnabled = !0, e.name = e.name.replace(o, "")), e && !s.test(e.name)) return i();
        e.name = e.name.replace(s, ""), y.resumeBootstrap = function(e) {
          return C(e, (function(e) {
            n.push(e)
          })), i()
        }, z(y.resumeDeferredBootstrap) && y.resumeDeferredBootstrap()
      }

      function Ee() {
        e.name = "NG_ENABLE_DEBUG_INFO!" + e.name, e.location.reload()
      }

      function Oe(e) {
        var t = y.element(e).injector();
        if (!t) throw b("test", "no injector found for element argument to getTestability");
        return t.get("$$testability")
      }
      var Ne = /[A-Z]/g;

      function Pe(e, t) {
        return t = t || "_", e.replace(Ne, (function(e, n) {
          return (n ? t : "") + e.toLowerCase()
        }))
      }
      var qe = !1;

      function Ie() {
        ft.legacyXHTMLReplacement = !0
      }

      function Le(e, t, n) {
        if (!e) throw b("areq", "Argument '{0}' is {1}", t || "?", n || "required");
        return e
      }

      function Re(e, t, n) {
        return n && G(e) && (e = e[e.length - 1]), Le(z(e), t, "not a function, got " + (e && "object" == typeof e ? e.constructor.name || "Object" :
          typeof e)), e
      }

      function je(e, t) {
        if ("hasOwnProperty" === e) throw b("badname", "hasOwnProperty is not a valid {0} name", t)
      }

      function _e(e) {
        for (var t, n = e[0], r = e[e.length - 1], i = 1; n !== r && (n = n.nextSibling); i++)(t || e[i] !== n) && (t || (t = a(h.call(e, 0, i))), t.push(
          n));
        return t || e
      }

      function Ve() {
        return Object.create(null)
      }

      function Ue(e) {
        if (null == e) return "";
        switch (typeof e) {
          case "string":
            break;
          case "number":
            e = "" + e;
            break;
          default:
            e = !R(e) || G(e) || B(e) ? ge(e) : e.toString()
        }
        return e
      }
      var He = 3;

      function Fe(e, t) {
        if (G(e)) {
          t = t || [];
          for (var n = 0, r = e.length; n < r; n++) t[n] = e[n]
        } else if (V(e))
          for (var i in t = t || {}, e) "$" === i.charAt(0) && "$" === i.charAt(1) || (t[i] = e[i]);
        return t || e
      }

      function Be(e, t) {
        return "function" == typeof e ? e.toString().replace(/ \{[\s\S]*$/, "") : j(e) ? "undefined" : "string" != typeof e ? function(e, t) {
          var n = [];
          return r(t) && (e = y.copy(e, null, t)), JSON.stringify(e, (function(e, t) {
            if (V(t = he(e, t))) {
              if (n.indexOf(t) >= 0) return "...";
              n.push(t)
            }
            return t
          }))
        }(e, t) : e
      }
      var Ge = {
        full: "1.8.2",
        major: 1,
        minor: 8,
        dot: 2,
        codeName: "meteoric-mining"
      };
      ft.expando = "ng339";
      var We = ft.cache = {},
        ze = 1;
      ft._data = function(e) {
        return this.cache[e[this.expando]] || {}
      };
      var Ye = /-([a-z])/g,
        Ke = /^-ms-/,
        Qe = {
          mouseleave: "mouseout",
          mouseenter: "mouseover"
        },
        Xe = i("jqLite");

      function Je(e, t) {
        return t.toUpperCase()
      }

      function Ze(e) {
        return e.replace(Ye, Je)
      }
      var et = /^<([\w-]+)\s*\/?>(?:<\/\1>|)$/,
        tt = /<|&#?\w+;/,
        nt = /<([\w:-]+)/,
        rt = /<(?!area|br|col|embed|hr|img|input|link|meta|param)(([\w:-]+)[^>]*)\/>/gi,
        it = {
          thead: ["table"],
          col: ["colgroup", "table"],
          tr: ["tbody", "table"],
          td: ["tr", "tbody", "table"]
        };
      it.tbody = it.tfoot = it.colgroup = it.caption = it.thead, it.th = it.td;
      var ot = {
        option: [1, '<select multiple="multiple">', "</select>"],
        _default: [0, "", ""]
      };
      for (var at in it) {
        var st = it[at],
          ut = st.slice().reverse();
        ot[at] = [ut.length, "<" + ut.join("><") + ">", "</" + st.join("></") + ">"]
      }

      function lt(e) {
        return !tt.test(e)
      }

      function ct(e) {
        var t = e.nodeType;
        return 1 === t || !t || 9 === t
      }

      function dt(t, n) {
        var r, i, a, s, u, l = n.createDocumentFragment(),
          c = [];
        if (lt(t)) c.push(n.createTextNode(t));
        else {
          if (r = l.appendChild(n.createElement("div")), i = (nt.exec(t) || ["", ""])[1].toLowerCase(), s = ft.legacyXHTMLReplacement ? t.replace(rt,
              "<$1></$2>") : t, o < 10)
            for (a = ot[i] || ot._default, r.innerHTML = a[1] + s + a[2], u = a[0]; u--;) r = r.firstChild;
          else {
            for (u = (a = it[i] || []).length; --u > -1;) r.appendChild(e.document.createElement(a[u])), r = r.firstChild;
            r.innerHTML = s
          }
          c = de(c, r.childNodes), (r = l.firstChild).textContent = ""
        }
        return l.textContent = "", l.innerHTML = "", C(c, (function(e) {
          l.appendChild(e)
        })), l
      }
      ot.optgroup = ot.option;
      var pt = e.Node.prototype.contains || function(e) {
        return !!(16 & this.compareDocumentPosition(e))
      };

      function ft(t) {
        if (t instanceof ft) return t;
        var n, r, i, o;
        if (H(t) && (t = ee(t), n = !0), !(this instanceof ft)) {
          if (n && "<" !== t.charAt(0)) throw Xe("nosel",
            "Looking up elements via selectors is not supported by jqLite! See: http://docs.angularjs.org/api/angular.element");
          return new ft(t)
        }
        n ? At(this, (r = t, i = i || e.document, (o = et.exec(r)) ? [i.createElement(o[1])] : (o = dt(r, i)) ? o.childNodes : [])) : z(t) ? Et(t) : At(
          this, t)
      }

      function ht(e) {
        return e.cloneNode(!0)
      }

      function gt(e, t) {
        !t && ct(e) && a.cleanData([e]), e.querySelectorAll && a.cleanData(e.querySelectorAll("*"))
      }

      function mt(e) {
        var t;
        for (t in e) return !1;
        return !0
      }

      function vt(e) {
        var t = e.ng339,
          n = t && We[t],
          r = n && n.events,
          i = n && n.data;
        i && !mt(i) || r && !mt(r) || (delete We[t], e.ng339 = void 0)
      }

      function $t(e, t, n, r) {
        if (_(r)) throw Xe("offargs", "jqLite#off() does not support the `selector` argument");
        var i = yt(e),
          o = i && i.events,
          a = i && i.handle;
        if (a) {
          if (t) {
            var s = function(t) {
              var r = o[t];
              _(n) && oe(r || [], n), _(n) && r && r.length > 0 || (e.removeEventListener(t, a), delete o[t])
            };
            C(t.split(" "), (function(e) {
              s(e), Qe[e] && s(Qe[e])
            }))
          } else
            for (t in o) "$destroy" !== t && e.removeEventListener(t, a), delete o[t];
          vt(e)
        }
      }

      function bt(e, t) {
        var n = e.ng339,
          r = n && We[n];
        r && (t ? delete r.data[t] : r.data = {}, vt(e))
      }

      function yt(e, t) {
        var n = e.ng339,
          r = n && We[n];
        return t && !r && (e.ng339 = n = ++ze, r = We[n] = {
          events: {},
          data: {},
          handle: void 0
        }), r
      }

      function wt(e, t, n) {
        if (ct(e)) {
          var r, i = _(n),
            o = !i && t && !V(t),
            a = !t,
            s = yt(e, !o),
            u = s && s.data;
          if (i) u[Ze(t)] = n;
          else {
            if (a) return u;
            if (o) return u && u[Ze(t)];
            for (r in t) u[Ze(r)] = t[r]
          }
        }
      }

      function xt(e, t) {
        return !!e.getAttribute && (" " + (e.getAttribute("class") || "") + " ").replace(/[\n\t]/g, " ").indexOf(" " + t + " ") > -1
      }

      function Ct(e, t) {
        if (t && e.setAttribute) {
          var n = (" " + (e.getAttribute("class") || "") + " ").replace(/[\n\t]/g, " "),
            r = n;
          C(t.split(" "), (function(e) {
            e = ee(e), r = r.replace(" " + e + " ", " ")
          })), r !== n && e.setAttribute("class", ee(r))
        }
      }

      function St(e, t) {
        if (t && e.setAttribute) {
          var n = (" " + (e.getAttribute("class") || "") + " ").replace(/[\n\t]/g, " "),
            r = n;
          C(t.split(" "), (function(e) {
            e = ee(e), -1 === r.indexOf(" " + e + " ") && (r += e + " ")
          })), r !== n && e.setAttribute("class", ee(r))
        }
      }

      function At(e, t) {
        if (t)
          if (t.nodeType) e[e.length++] = t;
          else {
            var n = t.length;
            if ("number" == typeof n && t.window !== t) {
              if (n)
                for (var r = 0; r < n; r++) e[e.length++] = t[r]
            } else e[e.length++] = t
          }
      }

      function Tt(e, t) {
        return kt(e, "$" + (t || "ngController") + "Controller")
      }

      function kt(e, t, n) {
        9 === e.nodeType && (e = e.documentElement);
        for (var r = G(t) ? t : [t]; e;) {
          for (var i = 0, o = r.length; i < o; i++)
            if (_(n = a.data(e, r[i]))) return n;
          e = e.parentNode || 11 === e.nodeType && e.host
        }
      }

      function Dt(e) {
        for (gt(e, !0); e.firstChild;) e.removeChild(e.firstChild)
      }

      function Mt(e, t) {
        t || gt(e);
        var n = e.parentNode;
        n && n.removeChild(e)
      }

      function Et(t) {
        function n() {
          e.document.removeEventListener("DOMContentLoaded", n), e.removeEventListener("load", n), t()
        }
        "complete" === e.document.readyState ? e.setTimeout(t) : (e.document.addEventListener("DOMContentLoaded", n), e.addEventListener("load", n))
      }
      var Ot = ft.prototype = {
          ready: Et,
          toString: function() {
            var e = [];
            return C(this, (function(t) {
              e.push("" + t)
            })), "[" + e.join(", ") + "]"
          },
          eq: function(e) {
            return a(e >= 0 ? this[e] : this[this.length + e])
          },
          length: 0,
          push: m,
          sort: [].sort,
          splice: [].splice
        },
        Nt = {};
      C("multiple,selected,checked,disabled,readOnly,required,open".split(","), (function(e) {
        Nt[p(e)] = e
      }));
      var Pt = {};
      C("input,select,option,textarea,button,form,details".split(","), (function(e) {
        Pt[e] = !0
      }));
      var qt = {
        ngMinlength: "minlength",
        ngMaxlength: "maxlength",
        ngMin: "min",
        ngMax: "max",
        ngPattern: "pattern",
        ngStep: "step"
      };

      function It(e, t) {
        var n = Nt[t.toLowerCase()];
        return n && Pt[re(e)] && n
      }

      function Lt(e, t, n) {
        n.call(e, t)
      }

      function Rt(e, t, n) {
        var r = t.relatedTarget;
        r && (r === e || pt.call(e, r)) || n.call(e, t)
      }

      function jt() {
        this.$get = function() {
          return M(ft, {
            hasClass: function(e, t) {
              return e.attr && (e = e[0]), xt(e, t)
            },
            addClass: function(e, t) {
              return e.attr && (e = e[0]), St(e, t)
            },
            removeClass: function(e, t) {
              return e.attr && (e = e[0]), Ct(e, t)
            }
          })
        }
      }

      function _t(e, t) {
        var n = e && e.$$hashKey;
        if (n) return "function" == typeof n && (n = e.$$hashKey()), n;
        var r = typeof e;
        return n = "function" === r || "object" === r && null !== e ? e.$$hashKey = r + ":" + (t || T)() : r + ":" + e
      }
      C({
        data: wt,
        removeData: bt,
        hasData: function(e) {
          for (var t in We[e.ng339]) return !0;
          return !1
        },
        cleanData: function(e) {
          for (var t = 0, n = e.length; t < n; t++) bt(e[t]), $t(e[t])
        }
      }, (function(e, t) {
        ft[t] = e
      })), C({
        data: wt,
        inheritedData: kt,
        scope: function(e) {
          return a.data(e, "$scope") || kt(e.parentNode || e, ["$isolateScope", "$scope"])
        },
        isolateScope: function(e) {
          return a.data(e, "$isolateScope") || a.data(e, "$isolateScopeNoTemplate")
        },
        controller: Tt,
        injector: function(e) {
          return kt(e, "$injector")
        },
        removeAttr: function(e, t) {
          e.removeAttribute(t)
        },
        hasClass: xt,
        css: function(e, t, n) {
          if (t = function(e) {
              return Ze(e.replace(Ke, "ms-"))
            }(t), !_(n)) return e.style[t];
          e.style[t] = n
        },
        attr: function(e, t, n) {
          var r, i = e.nodeType;
          if (i !== He && 2 !== i && 8 !== i && e.getAttribute) {
            var o = p(t),
              a = Nt[o];
            if (!_(n)) return r = e.getAttribute(t), a && null !== r && (r = o), null === r ? void 0 : r;
            null === n || !1 === n && a ? e.removeAttribute(t) : e.setAttribute(t, a ? o : n)
          }
        },
        prop: function(e, t, n) {
          if (!_(n)) return e[t];
          e[t] = n
        },
        text: function() {
          return e.$dv = "", e;

          function e(e, t) {
            if (j(t)) {
              var n = e.nodeType;
              return 1 === n || n === He ? e.textContent : ""
            }
            e.textContent = t
          }
        }(),
        val: function(e, t) {
          if (j(t)) {
            if (e.multiple && "select" === re(e)) {
              var n = [];
              return C(e.options, (function(e) {
                e.selected && n.push(e.value || e.text)
              })), n
            }
            return e.value
          }
          e.value = t
        },
        html: function(e, t) {
          if (j(t)) return e.innerHTML;
          gt(e, !0), e.innerHTML = t
        },
        empty: Dt
      }, (function(e, t) {
        ft.prototype[t] = function(t, n) {
          var r, i, o = this.length;
          if (e !== Dt && j(2 === e.length && e !== xt && e !== Tt ? t : n)) {
            if (V(t)) {
              for (r = 0; r < o; r++)
                if (e === wt) e(this[r], t);
                else
                  for (i in t) e(this[r], i, t[i]);
              return this
            }
            for (var a = e.$dv, s = j(a) ? Math.min(o, 1) : o, u = 0; u < s; u++) {
              var l = e(this[u], t, n);
              a = a ? a + l : l
            }
            return a
          }
          for (r = 0; r < o; r++) e(this[r], t, n);
          return this
        }
      })), C({
        removeData: bt,
        on: function(e, t, n, r) {
          if (_(r)) throw Xe("onargs", "jqLite#on() does not support the `selector` or `eventData` parameters");
          if (ct(e)) {
            var i = yt(e, !0),
              o = i.events,
              a = i.handle;
            a || (a = i.handle = function(e, t) {
              var n = function(n, r) {
                n.isDefaultPrevented = function() {
                  return n.defaultPrevented
                };
                var i = t[r || n.type],
                  o = i ? i.length : 0;
                if (o) {
                  if (j(n.immediatePropagationStopped)) {
                    var a = n.stopImmediatePropagation;
                    n.stopImmediatePropagation = function() {
                      n.immediatePropagationStopped = !0, n.stopPropagation && n.stopPropagation(), a && a.call(n)
                    }
                  }
                  n.isImmediatePropagationStopped = function() {
                    return !0 === n.immediatePropagationStopped
                  };
                  var s = i.specialHandlerWrapper || Lt;
                  o > 1 && (i = Fe(i));
                  for (var u = 0; u < o; u++) n.isImmediatePropagationStopped() || s(e, n, i[u])
                }
              };
              return n.elem = e, n
            }(e, o));
            for (var s = t.indexOf(" ") >= 0 ? t.split(" ") : [t], u = s.length, l = function(t, r, i) {
                var s = o[t];
                s || ((s = o[t] = []).specialHandlerWrapper = r, "$destroy" === t || i || e.addEventListener(t, a)), s.push(n)
              }; u--;) t = s[u], Qe[t] ? (l(Qe[t], Rt), l(t, void 0, !0)) : l(t)
          }
        },
        off: $t,
        one: function(e, t, n) {
          (e = a(e)).on(t, (function r() {
            e.off(t, n), e.off(t, r)
          })), e.on(t, n)
        },
        replaceWith: function(e, t) {
          var n, r = e.parentNode;
          gt(e), C(new ft(t), (function(t) {
            n ? r.insertBefore(t, n.nextSibling) : r.replaceChild(t, e), n = t
          }))
        },
        children: function(e) {
          var t = [];
          return C(e.childNodes, (function(e) {
            1 === e.nodeType && t.push(e)
          })), t
        },
        contents: function(e) {
          return e.contentDocument || e.childNodes || []
        },
        append: function(e, t) {
          var n = e.nodeType;
          if (1 === n || 11 === n)
            for (var r = 0, i = (t = new ft(t)).length; r < i; r++) {
              var o = t[r];
              e.appendChild(o)
            }
        },
        prepend: function(e, t) {
          if (1 === e.nodeType) {
            var n = e.firstChild;
            C(new ft(t), (function(t) {
              e.insertBefore(t, n)
            }))
          }
        },
        wrap: function(e, t) {
          var n, r, i;
          n = e, r = a(t).eq(0).clone()[0], (i = n.parentNode) && i.replaceChild(r, n), r.appendChild(n)
        },
        remove: Mt,
        detach: function(e) {
          Mt(e, !0)
        },
        after: function(e, t) {
          var n = e,
            r = e.parentNode;
          if (r)
            for (var i = 0, o = (t = new ft(t)).length; i < o; i++) {
              var a = t[i];
              r.insertBefore(a, n.nextSibling), n = a
            }
        },
        addClass: St,
        removeClass: Ct,
        toggleClass: function(e, t, n) {
          t && C(t.split(" "), (function(t) {
            var r = n;
            j(r) && (r = !xt(e, t)), (r ? St : Ct)(e, t)
          }))
        },
        parent: function(e) {
          var t = e.parentNode;
          return t && 11 !== t.nodeType ? t : null
        },
        next: function(e) {
          return e.nextElementSibling
        },
        find: function(e, t) {
          return e.getElementsByTagName ? e.getElementsByTagName(t) : []
        },
        clone: ht,
        triggerHandler: function(e, t, n) {
          var r, i, o, a = t.type || t,
            s = yt(e),
            u = s && s.events,
            l = u && u[a];
          l && (r = {
            preventDefault: function() {
              this.defaultPrevented = !0
            },
            isDefaultPrevented: function() {
              return !0 === this.defaultPrevented
            },
            stopImmediatePropagation: function() {
              this.immediatePropagationStopped = !0
            },
            isImmediatePropagationStopped: function() {
              return !0 === this.immediatePropagationStopped
            },
            stopPropagation: q,
            type: a,
            target: e
          }, t.type && (r = M(r, t)), i = Fe(l), o = n ? [r].concat(n) : [r], C(i, (function(t) {
            r.isImmediatePropagationStopped() || t.apply(e, o)
          })))
        }
      }, (function(e, t) {
        ft.prototype[t] = function(t, n, r) {
          for (var i, o = 0, s = this.length; o < s; o++) j(i) ? _(i = e(this[o], t, n, r)) && (i = a(i)) : At(i, e(this[o], t, n, r));
          return _(i) ? i : this
        }
      })), ft.prototype.bind = ft.prototype.on, ft.prototype.unbind = ft.prototype.off;
      var Vt = Object.create(null);

      function Ut() {
        this._keys = [], this._values = [], this._lastKey = NaN, this._lastIndex = -1
      }
      Ut.prototype = {
        _idx: function(e) {
          return e !== this._lastKey && (this._lastKey = e, this._lastIndex = this._keys.indexOf(e)), this._lastIndex
        },
        _transformKey: function(e) {
          return N(e) ? Vt : e
        },
        get: function(e) {
          e = this._transformKey(e);
          var t = this._idx(e);
          if (-1 !== t) return this._values[t]
        },
        has: function(e) {
          return e = this._transformKey(e), -1 !== this._idx(e)
        },
        set: function(e, t) {
          e = this._transformKey(e);
          var n = this._idx(e); - 1 === n && (n = this._lastIndex = this._keys.length), this._keys[n] = e, this._values[n] = t
        },
        delete: function(e) {
          e = this._transformKey(e);
          var t = this._idx(e);
          return -1 !== t && (this._keys.splice(t, 1), this._values.splice(t, 1), this._lastKey = NaN, this._lastIndex = -1, !0)
        }
      };
      var Ht = Ut,
        Ft = [function() {
          this.$get = [function() {
            return Ht
          }]
        }],
        Bt = /^([^(]+?)=>/,
        Gt = /^[^(]*\(\s*([^)]*)\)/m,
        Wt = /,/,
        zt = /^\s*(_?)(\S+?)\1\s*$/,
        Yt = /((\/\/.*$)|(\/\*[\s\S]*?\*\/))/gm,
        Kt = i("$injector");

      function Qt(e) {
        return Function.prototype.toString.call(e)
      }

      function Xt(e) {
        var t = Qt(e).replace(Yt, "");
        return t.match(Bt) || t.match(Gt)
      }

      function Jt(e, t) {
        t = !0 === t;
        var n = {},
          r = "Provider",
          i = [],
          a = new Ht,
          s = {
            $provide: {
              provider: h(g),
              factory: h(v),
              service: h((function(e, t) {
                return v(e, ["$injector", function(e) {
                  return e.instantiate(t)
                }])
              })),
              value: h((function(e, t) {
                return v(e, L(t), !1)
              })),
              constant: h((function(e, t) {
                je(e, "constant"), s[e] = t, c[e] = t
              })),
              decorator: function(e, t) {
                var n = l.get(e + r),
                  i = n.$get;
                n.$get = function() {
                  var e = p.invoke(i, n);
                  return p.invoke(t, null, {
                    $delegate: e
                  })
                }
              }
            }
          },
          l = s.$injector = b(s, (function(e, t) {
            throw y.isString(t) && i.push(t), Kt("unpr", "Unknown provider: {0}", i.join(" <- "))
          })),
          c = {},
          d = b(c, (function(e, t) {
            var n = l.get(e + r, t);
            return p.invoke(n.$get, n, void 0, e)
          })),
          p = d;
        s.$injectorProvider = {
          $get: L(d)
        }, p.modules = l.modules = Ve();
        var f = $(e);
        return (p = d.get("$injector")).strictDi = t, C(f, (function(e) {
          e && p.invoke(e)
        })), p.loadNewModules = function(e) {
          C($(e), (function(e) {
            e && p.invoke(e)
          }))
        }, p;

        function h(e) {
          return function(t, n) {
            if (!V(t)) return e(t, n);
            C(t, A(e))
          }
        }

        function g(e, t) {
          if (je(e, "service"), (z(t) || G(t)) && (t = l.instantiate(t)), !t.$get) throw Kt("pget", "Provider '{0}' must define $get factory method.", e);
          return s[e + r] = t
        }

        function m(e, t) {
          return function() {
            var n = p.invoke(t, this);
            if (j(n)) throw Kt("undef", "Provider '{0}' must return a value from $get factory method.", e);
            return n
          }
        }

        function v(e, t, n) {
          return g(e, {
            $get: !1 !== n ? m(e, t) : t
          })
        }

        function $(e) {
          Le(j(e) || G(e), "modulesToLoad", "not an array");
          var t, n = [];
          return C(e, (function(e) {
            if (!a.get(e)) {
              a.set(e, !0);
              try {
                H(e) ? (t = u(e), p.modules[e] = t, n = n.concat($(t.requires)).concat(t._runBlocks), r(t._invokeQueue), r(t._configBlocks)) : z(
                  e) || G(e) ? n.push(l.invoke(e)) : Re(e, "module")
              } catch (t) {
                throw G(e) && (e = e[e.length - 1]), t.message && t.stack && -1 === t.stack.indexOf(t.message) && (t = t.message + "\n" + t.stack),
                  Kt("modulerr", "Failed to instantiate module {0} due to:\n{1}", e, t.stack || t.message || t)
              }
            }

            function r(e) {
              var t, n;
              for (t = 0, n = e.length; t < n; t++) {
                var r = e[t],
                  i = l.get(r[0]);
                i[r[1]].apply(i, r[2])
              }
            }
          })), n
        }

        function b(e, a) {
          function u(t, r) {
            if (e.hasOwnProperty(t)) {
              if (e[t] === n) throw Kt("cdep", "Circular dependency found: {0}", t + " <- " + i.join(" <- "));
              return e[t]
            }
            try {
              return i.unshift(t), e[t] = n, e[t] = a(t, r), e[t]
            } catch (r) {
              throw e[t] === n && delete e[t], r
            } finally {
              i.shift()
            }
          }

          function l(e, n, r) {
            for (var i = [], o = Jt.$$annotate(e, t, r), a = 0, s = o.length; a < s; a++) {
              var l = o[a];
              if ("string" != typeof l) throw Kt("itkn", "Incorrect injection token! Expected service name as string, got {0}", l);
              i.push(n && n.hasOwnProperty(l) ? n[l] : u(l, r))
            }
            return i
          }
          return {
            invoke: function(e, t, n, r) {
              "string" == typeof n && (r = n, n = null);
              var i = l(e, n, r);
              return G(e) && (e = e[e.length - 1]),
                function(e) {
                  if (o || "function" != typeof e) return !1;
                  var t = e.$$ngIsClass;
                  return X(t) || (t = e.$$ngIsClass = /^class\b/.test(Qt(e))), t
                }(e) ? (i.unshift(null), new(Function.prototype.bind.apply(e, i))) : e.apply(t, i)
            },
            instantiate: function(e, t, n) {
              var r = G(e) ? e[e.length - 1] : e,
                i = l(e, t, n);
              return i.unshift(null), new(Function.prototype.bind.apply(r, i))
            },
            get: u,
            annotate: Jt.$$annotate,
            has: function(t) {
              return s.hasOwnProperty(t + r) || e.hasOwnProperty(t)
            }
          }
        }
      }

      function Zt() {
        var t = !0;
        this.disableAutoScrolling = function() {
          t = !1
        }, this.$get = ["$window", "$location", "$rootScope", function(n, r, i) {
          var o = n.document;

          function s(e) {
            if (e) {
              e.scrollIntoView();
              var t = function() {
                var e = u.yOffset;
                if (z(e)) e = e();
                else if (ne(e)) {
                  var t = e[0];
                  e = "fixed" !== n.getComputedStyle(t).position ? 0 : t.getBoundingClientRect().bottom
                } else F(e) || (e = 0);
                return e
              }();
              if (t) {
                var r = e.getBoundingClientRect().top;
                n.scrollBy(0, r - t)
              }
            } else n.scrollTo(0, 0)
          }

          function u(e) {
            var t, n, i;
            (e = H(e) ? e : F(e) ? e.toString() : r.hash()) ? (t = o.getElementById(e)) ? s(t): (n = o.getElementsByName(e), i = null, Array.prototype
              .some.call(n, (function(e) {
                if ("a" === re(e)) return i = e, !0
              })), (t = i) ? s(t) : "top" === e && s(null)): s(null)
          }
          return t && i.$watch((function() {
            return r.hash()
          }), (function(t, n) {
            var r, o;
            t === n && "" === t || (r = function() {
              i.$evalAsync(u)
            }, "complete" === (o = o || e).document.readyState ? o.setTimeout(r) : a(o).on("load", r))
          })), u
        }]
      }
      Jt.$$annotate = function(e, t, n) {
        var r, i;
        if ("function" == typeof e) {
          if (!(r = e.$inject)) {
            if (r = [], e.length) {
              if (t) throw H(n) && n || (n = e.name || function(e) {
                var t = Xt(e);
                return t ? "function(" + (t[1] || "").replace(/[\s\r\n]+/, " ") + ")" : "fn"
              }(e)), Kt("strictdi", "{0} is not using explicit annotation and cannot be invoked in strict mode", n);
              C(Xt(e)[1].split(Wt), (function(e) {
                e.replace(zt, (function(e, t, n) {
                  r.push(n)
                }))
              }))
            }
            e.$inject = r
          }
        } else G(e) ? (Re(e[i = e.length - 1], "fn"), r = e.slice(0, i)) : Re(e, "fn", !0);
        return r
      };
      var en = i("$animate"),
        tn = "ng-animate";

      function nn(e, t) {
        return e || t ? e ? t ? (G(e) && (e = e.join(" ")), G(t) && (t = t.join(" ")), e + " " + t) : e : t : ""
      }

      function rn(e) {
        return V(e) ? e : {}
      }
      var on = function() {
          this.$get = q
        },
        an = function() {
          var e = new Ht,
            t = [];
          this.$get = ["$$AnimateRunner", "$rootScope", function(n, r) {
            return {
              enabled: q,
              on: q,
              off: q,
              pin: q,
              push: function(a, s, u, l) {
                l && l(), (u = u || {}).from && a.css(u.from), u.to && a.css(u.to), (u.addClass || u.removeClass) && function(n, a, s) {
                  var u = e.get(n) || {},
                    l = i(u, a, !0),
                    c = i(u, s, !1);
                  (l || c) && (e.set(n, u), t.push(n), 1 === t.length && r.$$postDigest(o))
                }(a, u.addClass, u.removeClass);
                var c = new n;
                return c.complete(), c
              }
            };

            function i(e, t, n) {
              var r = !1;
              return t && C(t = H(t) ? t.split(" ") : G(t) ? t : [], (function(t) {
                t && (r = !0, e[t] = n)
              })), r
            }

            function o() {
              C(t, (function(t) {
                var n = e.get(t);
                if (n) {
                  var r = function(e) {
                      H(e) && (e = e.split(" "));
                      var t = Ve();
                      return C(e, (function(e) {
                        e.length && (t[e] = !0)
                      })), t
                    }(t.attr("class")),
                    i = "",
                    o = "";
                  C(n, (function(e, t) {
                    e !== !!r[t] && (e ? i += (i.length ? " " : "") + t : o += (o.length ? " " : "") + t)
                  })), C(t, (function(e) {
                    i && St(e, i), o && Ct(e, o)
                  })), e.delete(t)
                }
              })), t.length = 0
            }
          }]
        },
        sn = ["$provide", function(e) {
          var t = this,
            n = null,
            r = null;
          this.$$registeredAnimations = Object.create(null), this.register = function(n, r) {
            if (n && "." !== n.charAt(0)) throw en("notcsel", "Expecting class selector starting with '.' got '{0}'.", n);
            var i = n + "-animation";
            t.$$registeredAnimations[n.substr(1)] = i, e.factory(i, r)
          }, this.customFilter = function(e) {
            return 1 === arguments.length && (r = z(e) ? e : null), r
          }, this.classNameFilter = function(e) {
            if (1 === arguments.length && (n = e instanceof RegExp ? e : null)) {
              var t = new RegExp("[(\\s|\\/)]ng-animate[(\\s|\\/)]");
              if (t.test(n.toString())) throw n = null, en("nongcls",
                '$animateProvider.classNameFilter(regex) prohibits accepting a regex value which matches/contains the "{0}" CSS class.', tn)
            }
            return n
          }, this.$get = ["$$animateQueue", function(e) {
            function t(e, t, n) {
              if (n) {
                var r = function(e) {
                  for (var t = 0; t < e.length; t++) {
                    var n = e[t];
                    if (1 === n.nodeType) return n
                  }
                }(n);
                !r || r.parentNode || r.previousElementSibling || (n = null)
              }
              n ? n.after(e) : t.prepend(e)
            }
            return {
              on: e.on,
              off: e.off,
              pin: e.pin,
              enabled: e.enabled,
              cancel: function(e) {
                e.cancel && e.cancel()
              },
              enter: function(n, r, i, o) {
                return r = r && a(r), i = i && a(i), t(n, r = r || i.parent(), i), e.push(n, "enter", rn(o))
              },
              move: function(n, r, i, o) {
                return r = r && a(r), i = i && a(i), t(n, r = r || i.parent(), i), e.push(n, "move", rn(o))
              },
              leave: function(t, n) {
                return e.push(t, "leave", rn(n), (function() {
                  t.remove()
                }))
              },
              addClass: function(t, n, r) {
                return (r = rn(r)).addClass = nn(r.addclass, n), e.push(t, "addClass", r)
              },
              removeClass: function(t, n, r) {
                return (r = rn(r)).removeClass = nn(r.removeClass, n), e.push(t, "removeClass", r)
              },
              setClass: function(t, n, r, i) {
                return (i = rn(i)).addClass = nn(i.addClass, n), i.removeClass = nn(i.removeClass, r), e.push(t, "setClass", i)
              },
              animate: function(t, n, r, i, o) {
                return (o = rn(o)).from = o.from ? M(o.from, n) : n, o.to = o.to ? M(o.to, r) : r, i = i || "ng-inline-animate", o.tempClasses =
                  nn(o.tempClasses, i), e.push(t, "animate", o)
              }
            }
          }]
        }],
        un = function() {
          this.$get = ["$$rAF", function(e) {
            var t = [];

            function n(n) {
              t.push(n), t.length > 1 || e((function() {
                for (var e = 0; e < t.length; e++) t[e]();
                t = []
              }))
            }
            return function() {
              var e = !1;
              return n((function() {
                  e = !0
                })),
                function(t) {
                  e ? t() : n(t)
                }
            }
          }]
        },
        ln = function() {
          this.$get = ["$q", "$sniffer", "$$animateAsyncRun", "$$isDocumentHidden", "$timeout", function(e, t, n, r, i) {
            function o(e) {
              this.setHost(e);
              var t = n();
              this._doneCallbacks = [], this._tick = function(e) {
                r() ? function(e) {
                  i(e, 0, !1)
                }(e) : t(e)
              }, this._state = 0
            }
            return o.chain = function(e, t) {
              var n = 0;
              ! function r() {
                if (n === e.length) return void t(!0);
                e[n]((function(e) {
                  !1 !== e ? (n++, r()) : t(!1)
                }))
              }()
            }, o.all = function(e, t) {
              var n = 0,
                r = !0;

              function i(i) {
                r = r && i, ++n === e.length && t(r)
              }
              C(e, (function(e) {
                e.done(i)
              }))
            }, o.prototype = {
              setHost: function(e) {
                this.host = e || {}
              },
              done: function(e) {
                2 === this._state ? e() : this._doneCallbacks.push(e)
              },
              progress: q,
              getPromise: function() {
                if (!this.promise) {
                  var t = this;
                  this.promise = e((function(e, n) {
                    t.done((function(t) {
                      !1 === t ? n() : e()
                    }))
                  }))
                }
                return this.promise
              },
              then: function(e, t) {
                return this.getPromise().then(e, t)
              },
              catch: function(e) {
                return this.getPromise().catch(e)
              },
              finally: function(e) {
                return this.getPromise().finally(e)
              },
              pause: function() {
                this.host.pause && this.host.pause()
              },
              resume: function() {
                this.host.resume && this.host.resume()
              },
              end: function() {
                this.host.end && this.host.end(), this._resolve(!0)
              },
              cancel: function() {
                this.host.cancel && this.host.cancel(), this._resolve(!1)
              },
              complete: function(e) {
                var t = this;
                0 === t._state && (t._state = 1, t._tick((function() {
                  t._resolve(e)
                })))
              },
              _resolve: function(e) {
                2 !== this._state && (C(this._doneCallbacks, (function(t) {
                  t(e)
                })), this._doneCallbacks.length = 0, this._state = 2)
              }
            }, o
          }]
        },
        cn = function() {
          this.$get = ["$$rAF", "$q", "$$AnimateRunner", function(e, t, n) {
            return function(t, r) {
              var i = r || {};
              i.$$prepared || (i = ae(i)), i.cleanupStyles && (i.from = i.to = null), i.from && (t.css(i.from), i.from = null);
              var o, a = new n;
              return {
                start: s,
                end: s
              };

              function s() {
                return e((function() {
                  ! function() {
                    i.addClass && (t.addClass(i.addClass), i.addClass = null);
                    i.removeClass && (t.removeClass(i.removeClass), i.removeClass = null);
                    i.to && (t.css(i.to), i.to = null)
                  }(), o || a.complete(), o = !0
                })), a
              }
            }
          }]
        };

      function dn(e, t, n, r, i) {
        var o = this,
          s = e.location,
          u = e.history,
          l = e.setTimeout,
          c = e.clearTimeout,
          d = {},
          p = i(n);
        o.isMock = !1, o.$$completeOutstandingRequest = p.completeTask, o.$$incOutstandingRequestCount = p.incTaskCount, o
          .notifyWhenNoOutstandingRequests = p.notifyWhenNoPendingTasks;
        var f, h, g = s.href,
          m = t.find("base"),
          v = null,
          $ = r.history ? function() {
            try {
              return u.state
            } catch (e) {}
          } : q;
        S(), o.url = function(t, n, i) {
          if (j(i) && (i = null), s !== e.location && (s = e.location), u !== e.history && (u = e.history), t) {
            var a = h === i;
            if (t = di(t).href, g === t && (!r.history || a)) return o;
            var l = g && cr(g) === cr(t);
            return g = t, h = i, !r.history || l && a ? (l || (v = t), n ? s.replace(t) : l ? s.hash = function(e) {
              var t = e.indexOf("#");
              return -1 === t ? "" : e.substr(t)
            }(t) : s.href = t, s.href !== t && (v = t)) : (u[n ? "replaceState" : "pushState"](i, "", t), S()), v && (v = t), o
          }
          return function(e) {
            return e.replace(/#$/, "")
          }(v || s.href)
        }, o.state = function() {
          return f
        };
        var b = [],
          y = !1;

        function w() {
          v = null, A()
        }
        var x = null;

        function S() {
          ue(f = j(f = $()) ? null : f, x) && (f = x), x = f, h = f
        }

        function A() {
          var e = h;
          S(), g === o.url() && e === f || (g = o.url(), h = f, C(b, (function(e) {
            e(o.url(), f)
          })))
        }
        o.onUrlChange = function(t) {
          return y || (r.history && a(e).on("popstate", w), a(e).on("hashchange", w), y = !0), b.push(t), t
        }, o.$$applicationDestroyed = function() {
          a(e).off("hashchange popstate", w)
        }, o.$$checkUrlChange = A, o.baseHref = function() {
          var e = m.attr("href");
          return e ? e.replace(/^(https?:)?\/\/[^/]*/, "") : ""
        }, o.defer = function(e, t, n) {
          var r;
          return t = t || 0, n = n || p.DEFAULT_TASK_TYPE, p.incTaskCount(n), r = l((function() {
            delete d[r], p.completeTask(e, n)
          }), t), d[r] = n, r
        }, o.defer.cancel = function(e) {
          if (d.hasOwnProperty(e)) {
            var t = d[e];
            return delete d[e], c(e), p.completeTask(q, t), !0
          }
          return !1
        }
      }

      function pn() {
        this.$get = ["$window", "$log", "$sniffer", "$document", "$$taskTrackerFactory", function(e, t, n, r, i) {
          return new dn(e, r, t, n, i)
        }]
      }

      function fn() {
        this.$get = function() {
          var e = {};

          function t(t, n) {
            if (t in e) throw i("$cacheFactory")("iid", "CacheId '{0}' is already taken!", t);
            var r = 0,
              o = M({}, n, {
                id: t
              }),
              a = Ve(),
              s = n && n.capacity || Number.MAX_VALUE,
              u = Ve(),
              l = null,
              c = null;
            return e[t] = {
              put: function(e, t) {
                if (!j(t)) {
                  if (s < Number.MAX_VALUE) d(u[e] || (u[e] = {
                    key: e
                  }));
                  return e in a || r++, a[e] = t, r > s && this.remove(c.key), t
                }
              },
              get: function(e) {
                if (s < Number.MAX_VALUE) {
                  var t = u[e];
                  if (!t) return;
                  d(t)
                }
                return a[e]
              },
              remove: function(e) {
                if (s < Number.MAX_VALUE) {
                  var t = u[e];
                  if (!t) return;
                  t === l && (l = t.p), t === c && (c = t.n), p(t.n, t.p), delete u[e]
                }
                e in a && (delete a[e], r--)
              },
              removeAll: function() {
                a = Ve(), r = 0, u = Ve(), l = c = null
              },
              destroy: function() {
                a = null, o = null, u = null, delete e[t]
              },
              info: function() {
                return M({}, o, {
                  size: r
                })
              }
            };

            function d(e) {
              e !== l && (c ? c === e && (c = e.n) : c = e, p(e.n, e.p), p(e, l), (l = e).n = null)
            }

            function p(e, t) {
              e !== t && (e && (e.p = t), t && (t.n = e))
            }
          }
          return t.info = function() {
            var t = {};
            return C(e, (function(e, n) {
              t[n] = e.info()
            })), t
          }, t.get = function(t) {
            return e[t]
          }, t
        }
      }

      function hn() {
        this.$get = ["$cacheFactory", function(e) {
          return e("templates")
        }]
      }
      var gn = i("$compile");
      var mn = new function() {};

      function vn(t, n) {
        var r = {},
          i = "Directive",
          s = /^\s*directive:\s*([\w-]+)\s+(.*)$/,
          u = /(([\w-]+)(?::([^;]+))?;?)/,
          l = function(e) {
            var t, n = {},
              r = e.split(",");
            for (t = 0; t < r.length; t++) n[r[t]] = !0;
            return n
          }("ngSrc,ngSrcset,src,srcset"),
          c = /^(?:(\^\^?)?(\?)?(\^\^?)?)?/,
          f = /^(on[a-z]+|formaction)$/,
          h = Ve();

        function g(e, t, n) {
          var r = /^([@&]|[=<](\*?))(\??)\s*([\w$]*)$/,
            i = Ve();
          return C(e, (function(e, o) {
            if ((e = e.trim()) in h) i[o] = h[e];
            else {
              var a = e.match(r);
              if (!a) throw gn("iscp", "Invalid {3} for directive '{0}'. Definition: {... {1}: '{2}' ...}", t, o, e, n ?
                "controller bindings definition" : "isolate scope definition");
              i[o] = {
                mode: a[1][0],
                collection: "*" === a[2],
                optional: "?" === a[3],
                attrName: a[4] || o
              }, a[4] && (h[e] = i[o])
            }
          })), i
        }

        function m(e, t) {
          var n = {
            isolateScope: null,
            bindToController: null
          };
          if (V(e.scope) && (!0 === e.bindToController ? (n.bindToController = g(e.scope, t, !0), n.isolateScope = {}) : n.isolateScope = g(e.scope, t, !
              1)), V(e.bindToController) && (n.bindToController = g(e.bindToController, t, !0)), n.bindToController && !e.controller) throw gn("noctrl",
            "Cannot bind to controller without directive '{0}'s controller.", t);
          return n
        }
        this.directive = function e(n, o) {
          return Le(n, "name"), je(n, "directive"), H(n) ? (! function(e) {
            var t = e.charAt(0);
            if (!t || t !== p(t)) throw gn("baddir", "Directive/Component name '{0}' is invalid. The first character must be a lowercase letter",
            e);
            if (e !== e.trim()) throw gn("baddir",
              "Directive/Component name '{0}' is invalid. The name should not contain leading or trailing whitespaces", e)
          }(n), Le(o, "directiveFactory"), r.hasOwnProperty(n) || (r[n] = [], t.factory(n + i, ["$injector", "$exceptionHandler", function(e, t) {
            var i = [];
            return C(r[n], (function(r, o) {
              try {
                var a = e.invoke(r);
                z(a) ? a = {
                    compile: L(a)
                  } : !a.compile && a.link && (a.compile = L(a.link)), a.priority = a.priority || 0, a.index = o, a.name = a.name || n,
                  a.require = function(e) {
                    var t = e.require || e.controller && e.name;
                    return !G(t) && V(t) && C(t, (function(e, n) {
                      var r = e.match(c);
                      e.substring(r[0].length) || (t[n] = r[0] + n)
                    })), t
                  }(a), a.restrict = function(e, t) {
                    if (e && (!H(e) || !/[EACM]/.test(e))) throw gn("badrestrict",
                      "Restrict property '{0}' of directive '{1}' is invalid", e, t);
                    return e || "EA"
                  }(a.restrict, n), a.$$moduleName = r.$$moduleName, i.push(a)
              } catch (e) {
                t(e)
              }
            })), i
          }])), r[n].push(o)) : C(n, A(e)), this
        }, this.component = function e(t, n) {
          if (!H(t)) return C(t, A(fe(this, e))), this;
          var r = n.controller || function() {};

          function i(e) {
            function t(t) {
              return z(t) || G(t) ? function(n, r) {
                return e.invoke(t, this, {
                  $element: n,
                  $attrs: r
                })
              } : t
            }
            var i = n.template || n.templateUrl ? n.template : "",
              o = {
                controller: r,
                controllerAs: Tn(n.controller) || n.controllerAs || "$ctrl",
                template: t(i),
                templateUrl: t(n.templateUrl),
                transclude: n.transclude,
                scope: {},
                bindToController: n.bindings || {},
                restrict: "E",
                require: n.require
              };
            return C(n, (function(e, t) {
              "$" === t.charAt(0) && (o[t] = e)
            })), o
          }
          return C(n, (function(e, t) {
            "$" === t.charAt(0) && (i[t] = e, z(r) && (r[t] = e))
          })), i.$inject = ["$injector"], this.directive(t, i)
        }, this.aHrefSanitizationTrustedUrlList = function(e) {
          return _(e) ? (n.aHrefSanitizationTrustedUrlList(e), this) : n.aHrefSanitizationTrustedUrlList()
        }, Object.defineProperty(this, "aHrefSanitizationWhitelist", {
          get: function() {
            return this.aHrefSanitizationTrustedUrlList
          },
          set: function(e) {
            this.aHrefSanitizationTrustedUrlList = e
          }
        }), this.imgSrcSanitizationTrustedUrlList = function(e) {
          return _(e) ? (n.imgSrcSanitizationTrustedUrlList(e), this) : n.imgSrcSanitizationTrustedUrlList()
        }, Object.defineProperty(this, "imgSrcSanitizationWhitelist", {
          get: function() {
            return this.imgSrcSanitizationTrustedUrlList
          },
          set: function(e) {
            this.imgSrcSanitizationTrustedUrlList = e
          }
        });
        var $ = !0;
        this.debugInfoEnabled = function(e) {
          return _(e) ? ($ = e, this) : $
        };
        var b = !1;
        this.strictComponentBindingsEnabled = function(e) {
          return _(e) ? (b = e, this) : b
        };
        var y = 10;
        this.onChangesTtl = function(e) {
          return arguments.length ? (y = e, this) : y
        };
        var w = !0;
        this.commentDirectivesEnabled = function(e) {
          return arguments.length ? (w = e, this) : w
        };
        var x = !0;
        this.cssClassDirectivesEnabled = function(e) {
          return arguments.length ? (x = e, this) : x
        };
        var S = Ve();
        this.addPropertySecurityContext = function(e, t, n) {
            var r = e.toLowerCase() + "|" + t.toLowerCase();
            if (r in S && S[r] !== n) throw gn("ctxoverride", "Property context '{0}.{1}' already set to '{2}', cannot override to '{3}'.", e, t, S[r],
            n);
            return S[r] = n, this
          },
          function() {
            function e(e, t) {
              C(t, (function(t) {
                S[t.toLowerCase()] = e
              }))
            }
            e(zr.HTML, ["iframe|srcdoc", "*|innerHTML", "*|outerHTML"]), e(zr.CSS, ["*|style"]), e(zr.URL, ["area|href", "area|ping", "a|href", "a|ping",
              "blockquote|cite", "body|background", "del|cite", "input|src", "ins|cite", "q|cite"
            ]), e(zr.MEDIA_URL, ["audio|src", "img|src", "img|srcset", "source|src", "source|srcset", "track|src", "video|src", "video|poster"]), e(zr
              .RESOURCE_URL, ["*|formAction", "applet|code", "applet|codebase", "base|href", "embed|src", "frame|src", "form|action", "head|profile",
                "html|manifest", "iframe|src", "link|href", "media|src", "object|codebase", "object|data", "script|src"
              ])
          }(), this.$get = ["$injector", "$interpolate", "$exceptionHandler", "$templateRequest", "$parse", "$controller", "$rootScope", "$sce",
            "$animate",
            function(t, n, h, g, A, T, k, D, E) {
              var O, N = /^\w/,
                L = e.document.createElement("div"),
                R = w,
                _ = x,
                U = y;

              function F() {
                try {
                  if (!--U) throw O = void 0, gn("infchng", "{0} $onChanges() iterations reached. Aborting!\n", y);
                  k.$apply((function() {
                    for (var e = 0, t = O.length; e < t; ++e) try {
                      O[e]()
                    } catch (e) {
                      h(e)
                    }
                    O = void 0
                  }))
                } finally {
                  U++
                }
              }

              function B(e, t) {
                if (!e) return e;
                if (!H(e)) throw gn("srcset", 'Can\'t pass trusted values to `{0}`: "{1}"', t, e.toString());
                for (var n = "", r = ee(e), i = /\s/.test(r) ? /(\s+\d+x\s*,|\s+\d+w\s*,|\s+,|,\s+)/ : /(,)/, o = r.split(i), a = Math.floor(o.length /
                    2), s = 0; s < a; s++) {
                  var u = 2 * s;
                  n += D.getTrustedMediaUrl(ee(o[u])), n += " " + ee(o[u + 1])
                }
                var l = ee(o[2 * s]).split(/\s/);
                return n += D.getTrustedMediaUrl(ee(l[0])), 2 === l.length && (n += " " + ee(l[1])), n
              }

              function Y(e, t) {
                if (t) {
                  var n, r, i, o = Object.keys(t);
                  for (n = 0, r = o.length; n < r; n++) this[i = o[n]] = t[i]
                } else this.$attr = {};
                this.$$element = e
              }

              function K(e, t) {
                try {
                  e.addClass(t)
                } catch (e) {}
              }
              Y.prototype = {
                $normalize: wn,
                $addClass: function(e) {
                  e && e.length > 0 && E.addClass(this.$$element, e)
                },
                $removeClass: function(e) {
                  e && e.length > 0 && E.removeClass(this.$$element, e)
                },
                $updateClass: function(e, t) {
                  var n = xn(e, t);
                  n && n.length && E.addClass(this.$$element, n);
                  var r = xn(t, e);
                  r && r.length && E.removeClass(this.$$element, r)
                },
                $set: function(e, t, n, r) {
                  var i = It(this.$$element[0], e),
                    o = qt[e],
                    a = e;
                  i ? (this.$$element.prop(e, t), r = i) : o && (this[o] = t, a = o), this[e] = t, r ? this.$attr[e] = r : (r = this.$attr[e]) || (
                      this.$attr[e] = r = Pe(e, "-")), "img" === re(this.$$element) && "srcset" === e && (this[e] = t = B(t,
                    "$set('srcset', value)")), !1 !== n && (null === t || j(t) ? this.$$element.removeAttr(r) : N.test(r) ? i && !1 === t ? this
                      .$$element.removeAttr(r) : this.$$element.attr(r, t) : function(e, t, n) {
                        L.innerHTML = "<span " + t + ">";
                        var r = L.firstChild.attributes,
                          i = r[0];
                        r.removeNamedItem(i.name), i.value = n, e.attributes.setNamedItem(i)
                      }(this.$$element[0], r, t));
                  var s = this.$$observers;
                  s && C(s[a], (function(e) {
                    try {
                      e(t)
                    } catch (e) {
                      h(e)
                    }
                  }))
                },
                $observe: function(e, t) {
                  var n = this,
                    r = n.$$observers || (n.$$observers = Ve()),
                    i = r[e] || (r[e] = []);
                  return i.push(t), k.$evalAsync((function() {
                      i.$$inter || !n.hasOwnProperty(e) || j(n[e]) || t(n[e])
                    })),
                    function() {
                      oe(i, t)
                    }
                }
              };
              var J = n.startSymbol(),
                Z = n.endSymbol(),
                te = "{{" === J && "}}" === Z ? I : function(e) {
                  return e.replace(/\{\{/g, J).replace(/}}/g, Z)
                },
                ne = /^ng(Attr|Prop|On)([A-Z].*)$/,
                ie = /^(.+)Start$/;
              return ae.$$addBindingInfo = $ ? function(e, t) {
                var n = e.data("$binding") || [];
                G(t) ? n = n.concat(t) : n.push(t), e.data("$binding", n)
              } : q, ae.$$addBindingClass = $ ? function(e) {
                K(e, "ng-binding")
              } : q, ae.$$addScopeInfo = $ ? function(e, t, n, r) {
                var i = n ? r ? "$isolateScopeNoTemplate" : "$isolateScope" : "$scope";
                e.data(i, t)
              } : q, ae.$$addScopeClass = $ ? function(e, t) {
                K(e, t ? "ng-isolate-scope" : "ng-scope")
              } : q, ae.$$createComment = function(t, n) {
                var r = "";
                return $ && (r = " " + (t || "") + ": ", n && (r += n + " ")), e.document.createComment(r)
              }, ae;

              function ae(e, t, n, r, i) {
                e instanceof a || (e = a(e));
                var o = le(e, t, e, n, r, i);
                ae.$$addScopeClass(e);
                var s = null;
                return function(t, n, r) {
                  if (!e) throw gn("multilink", "This element has already been linked.");
                  Le(t, "scope"), i && i.needsNewScope && (t = t.$parent.$new());
                  var u, l, c, d = (r = r || {}).parentBoundTranscludeFn,
                    p = r.transcludeControllers,
                    f = r.futureParentElement;
                  if (d && d.$$boundTransclude && (d = d.$$boundTransclude), s || (l = (u = f) && u[0], s = l && "foreignobject" !== re(l) && v.call(l)
                      .match(/SVG/) ? "svg" : "html"), c = "html" !== s ? a(De(s, a("<div></div>").append(e).html())) : n ? Ot.clone.call(e) : e, p)
                    for (var h in p) c.data("$" + h + "Controller", p[h].instance);
                  return ae.$$addScopeInfo(c, t), n && n(c, t), o && o(t, c, c, d), n || (e = o = null), c
                }
              }

              function le(e, t, n, r, i, s) {
                for (var u, l, c, d, p, f, h, g = [], m = G(e) || e instanceof a, v = 0; v < e.length; v++) u = new Y, 11 === o && ce(e, v, m), (c = (l =
                  he(e[v], [], u, 0 === v ? r : void 0, i)).length ? $e(l, e[v], u, t, n, null, [], [], s) : null) && c.scope && ae.$$addScopeClass(u
                  .$$element), p = c && c.terminal || !(d = e[v].childNodes) || !d.length ? null : le(d, c ? (c.transcludeOnThisElement || !c
                  .templateOnThisElement) && c.transclude : t), (c || p) && (g.push(v, c, p), f = !0, h = h || c), s = null;
                return f ? function(e, n, r, i) {
                  var o, s, u, l, c, d, p, f, m;
                  if (h) {
                    var v = n.length;
                    for (m = new Array(v), c = 0; c < g.length; c += 3) p = g[c], m[p] = n[p]
                  } else m = n;
                  for (c = 0, d = g.length; c < d;) u = m[g[c++]], o = g[c++], s = g[c++], o ? (o.scope ? (l = e.$new(), ae.$$addScopeInfo(a(u), l)) :
                    l = e, f = o.transcludeOnThisElement ? de(e, o.transclude, i) : !o.templateOnThisElement && i ? i : !i && t ? de(e, t) : null, o(
                      s, l, u, r, f)) : s && s(e, u.childNodes, void 0, i)
                } : null
              }

              function ce(e, t, n) {
                var r, i = e[t],
                  o = i.parentNode;
                if (i.nodeType === He)
                  for (;
                    (r = o ? i.nextSibling : e[t + 1]) && r.nodeType === He;) i.nodeValue = i.nodeValue + r.nodeValue, r.parentNode && r.parentNode
                    .removeChild(r), n && r === e[t + 1] && e.splice(t + 1, 1)
              }

              function de(e, t, n) {
                function r(r, i, o, a, s) {
                  return r || ((r = e.$new(!1, s)).$$transcluded = !0), t(r, i, {
                    parentBoundTranscludeFn: n,
                    transcludeControllers: o,
                    futureParentElement: a
                  })
                }
                var i = r.$$slots = Ve();
                for (var o in t.$$slots) t.$$slots[o] ? i[o] = de(e, t.$$slots[o], n) : i[o] = null;
                return r
              }

              function he(e, t, r, i, o) {
                var a, l, c, d = e.nodeType,
                  p = r.$attr;
                switch (d) {
                  case 1:
                    xe(t, wn(l = re(e)), "E", i, o);
                    for (var f, h, g, m, v, $ = e.attributes, b = 0, y = $ && $.length; b < y; b++) {
                      var w, x = !1,
                        C = !1,
                        S = !1,
                        A = !1,
                        T = !1;
                      h = (f = $[b]).name, m = f.value, (v = (g = wn(h.toLowerCase())).match(ne)) ? (S = "Attr" === v[1], A = "Prop" === v[1], T =
                          "On" === v[1], h = h.replace(bn, "").toLowerCase().substr(4 + v[1].length).replace(/_(.)/g, (function(e, t) {
                            return t.toUpperCase()
                          }))) : (w = g.match(ie)) && Ce(w[1]) && (x = h, C = h.substr(0, h.length - 5) + "end", h = h.substr(0, h.length - 6)), A || T ?
                        (r[g] = m, p[g] = f.name, A ? Ee(e, t, g, h) : Oe(t, g, h)) : (p[g = wn(h.toLowerCase())] = h, !S && r.hasOwnProperty(g) || (r[
                          g] = m, It(e, g) && (r[g] = !0)), Ne(e, t, m, g, S), xe(t, g, "A", i, o, x, C))
                    }
                    if ("input" === l && "hidden" === e.getAttribute("type") && e.setAttribute("autocomplete", "off"), !_) break;
                    if (V(c = e.className) && (c = c.animVal), H(c) && "" !== c)
                      for (; a = u.exec(c);) xe(t, g = wn(a[2]), "C", i, o) && (r[g] = ee(a[3])), c = c.substr(a.index + a[0].length);
                    break;
                  case He:
                    ! function(e, t) {
                      var r = n(t, !0);
                      r && e.push({
                        priority: 0,
                        compile: function(e) {
                          var t = e.parent(),
                            n = !!t.length;
                          return n && ae.$$addBindingClass(t),
                            function(e, t) {
                              var i = t.parent();
                              n || ae.$$addBindingClass(i), ae.$$addBindingInfo(i, r.expressions), e.$watch(r, (function(e) {
                                t[0].nodeValue = e
                              }))
                            }
                        }
                      })
                    }(t, e.nodeValue);
                    break;
                  case 8:
                    if (!R) break;
                    ! function(e, t, n, r, i) {
                      try {
                        var o = s.exec(e.nodeValue);
                        if (o) {
                          var a = wn(o[1]);
                          xe(t, a, "M", r, i) && (n[a] = ee(o[2]))
                        }
                      } catch (e) {}
                    }(e, t, r, i, o)
                }
                return t.sort(Te), t
              }

              function ge(e, t, n) {
                var r = [],
                  i = 0;
                if (t && e.hasAttribute && e.hasAttribute(t))
                  do {
                    if (!e) throw gn("uterdir", "Unterminated attribute, found '{0}' but no matching '{1}' found.", t, n);
                    1 === e.nodeType && (e.hasAttribute(t) && i++, e.hasAttribute(n) && i--), r.push(e), e = e.nextSibling
                  } while (i > 0);
                else r.push(e);
                return a(r)
              }

              function me(e, t, n) {
                return function(r, i, o, a, s) {
                  return i = ge(i[0], t, n), e(r, i, o, a, s)
                }
              }

              function ve(e, t, n, r, i, o) {
                var a;
                return e ? ae(t, n, r, i, o) : function() {
                  return a || (a = ae(t, n, r, i, o), t = n = o = null), a.apply(this, arguments)
                }
              }

              function $e(t, n, r, i, o, s, u, l, c) {
                c = c || {};
                for (var d, p, f, g, m, v = -Number.MAX_VALUE, $ = c.newScopeDirective, b = c.controllerDirectives, y = c.newIsolateScopeDirective, w = c
                    .templateDirective, x = c.nonTlbTranscludeDirective, S = !1, A = !1, k = c.hasElementTranscludeDirective, D = r.$$element = a(n), E =
                    s, O = i, N = !1, P = !1, q = 0, I = t.length; q < I; q++) {
                  var L = (d = t[q]).$$start,
                    R = d.$$end;
                  if (L && (D = ge(n, L, R)), f = void 0, v > d.priority) break;
                  if ((m = d.scope) && (d.templateUrl || (V(m) ? (ke("new/isolated scope", y || $, d, D), y = d) : ke("new/isolated scope", y, d, D)), $ =
                      $ || d), p = d.name, !N && (d.replace && (d.templateUrl || d.template) || d.transclude && !d.$$tlb)) {
                    for (var _, U = q + 1; _ = t[U++];)
                      if (_.transclude && !_.$$tlb || _.replace && (_.templateUrl || _.template)) {
                        P = !0;
                        break
                      } N = !0
                  }
                  if (!d.templateUrl && d.controller && (b = b || Ve(), ke("'" + p + "' controller", b[p], d, D), b[p] = d), m = d.transclude)
                    if (S = !0, d.$$tlb || (ke("transclusion", x, d, D), x = d), "element" === m) k = !0, v = d.priority, f = D, D = r.$$element = a(ae
                      .$$createComment(p, r[p])), n = D[0], qe(o, pe(f), n), O = ve(P, f, i, v, E && E.name, {
                      nonTlbTranscludeDirective: x
                    });
                    else {
                      var H = Ve();
                      if (V(m)) {
                        f = e.document.createDocumentFragment();
                        var F = Ve(),
                          B = Ve();
                        for (var W in C(m, (function(e, t) {
                            var n = "?" === e.charAt(0);
                            e = n ? e.substring(1) : e, F[e] = t, H[t] = null, B[t] = n
                          })), C(D.contents(), (function(t) {
                            var n = F[wn(re(t))];
                            n ? (B[n] = !0, H[n] = H[n] || e.document.createDocumentFragment(), H[n].appendChild(t)) : f.appendChild(t)
                          })), C(B, (function(e, t) {
                            if (!e) throw gn("reqslot", "Required transclusion slot `{0}` was not filled.", t)
                          })), H)
                          if (H[W]) {
                            var K = a(H[W].childNodes);
                            H[W] = ve(P, K, i)
                          } f = a(f.childNodes)
                      } else f = a(ht(n)).contents();
                      D.empty(), (O = ve(P, f, i, void 0, void 0, {
                        needsNewScope: d.$$isolateScope || d.$$newScope
                      })).$$slots = H
                    } if (d.template)
                    if (A = !0, ke("template", w, d, D), w = d, m = z(d.template) ? d.template(D, r) : d.template, m = te(m), d.replace) {
                      if (E = d, f = lt(m) ? [] : Cn(De(d.templateNamespace, ee(m))), n = f[0], 1 !== f.length || 1 !== n.nodeType) throw gn("tplrt",
                        "Template for directive '{0}' must have exactly one root element. {1}", p, "");
                      qe(o, D, n);
                      var X = {
                          $attr: {}
                        },
                        J = he(n, [], X),
                        Z = t.splice(q + 1, t.length - (q + 1));
                      (y || $) && ye(J, y, $), t = t.concat(J).concat(Z), Se(r, X), I = t.length
                    } else D.html(m);
                  if (d.templateUrl) A = !0, ke("template", w, d, D), w = d, d.replace && (E = d), oe = Ae(t.splice(q, t.length - q), D, r, o, S && O, u,
                    l, {
                      controllerDirectives: b,
                      newScopeDirective: $ !== d && $,
                      newIsolateScopeDirective: y,
                      templateDirective: w,
                      nonTlbTranscludeDirective: x
                    }), I = t.length;
                  else if (d.compile) try {
                    g = d.compile(D, r, O);
                    var ne = d.$$originalDirective || d;
                    z(g) ? ie(null, fe(ne, g), L, R) : g && ie(fe(ne, g.pre), fe(ne, g.post), L, R)
                  } catch (e) {
                    h(e, we(D))
                  }
                  d.terminal && (oe.terminal = !0, v = Math.max(v, d.priority))
                }
                return oe.scope = $ && !0 === $.scope, oe.transcludeOnThisElement = S, oe.templateOnThisElement = A, oe.transclude = O, c
                  .hasElementTranscludeDirective = k, oe;

                function ie(e, t, n, r) {
                  e && (n && (e = me(e, n, r)), e.require = d.require, e.directiveName = p, (y === d || d.$$isolateScope) && (e = Ie(e, {
                    isolateScope: !0
                  })), u.push(e)), t && (n && (t = me(t, n, r)), t.require = d.require, t.directiveName = p, (y === d || d.$$isolateScope) && (t = Ie(
                  t, {
                    isolateScope: !0
                  })), l.push(t))
                }

                function oe(e, t, i, o, s) {
                  var c, d, p, f, g, m, v, x, S, A;
                  for (var D in n === i ? (S = r, x = r.$$element) : S = new Y(x = a(i), r), g = t, y ? f = t.$new(!0) : $ && (g = t.$parent), s && ((v =
                      function(e, t, n, r) {
                        var i;
                        Q(e) || (r = n, n = t, t = e, e = void 0);
                        k && (i = m);
                        n || (n = k ? x.parent() : x);
                        if (!r) return s(e, t, i, n, P);
                        var o = s.$$slots[r];
                        if (o) return o(e, t, i, n, P);
                        if (j(o)) throw gn("noslot", 'No parent directive that requires a transclusion with slot name "{0}". Element: {1}', r, we(x))
                      }).$$boundTransclude = s, v.isSlotFilled = function(e) {
                      return !!s.$$slots[e]
                    }), b && (m = function(e, t, n, r, i, o, a) {
                      var s = Ve();
                      for (var u in r) {
                        var l = r[u],
                          c = {
                            $scope: l === a || l.$$isolateScope ? i : o,
                            $element: e,
                            $attrs: t,
                            $transclude: n
                          },
                          d = l.controller;
                        "@" === d && (d = t[l.name]);
                        var p = T(d, c, !0, l.controllerAs);
                        s[l.name] = p, e.data("$" + l.name + "Controller", p.instance)
                      }
                      return s
                    }(x, S, v, b, f, t, y)), y && (ae.$$addScopeInfo(x, f, !0, !(w && (w === y || w === y.$$originalDirective))), ae.$$addScopeClass(x,
                      !0), f.$$isolateBindings = y.$$isolateBindings, (A = _e(t, S, f, f.$$isolateBindings, y)).removeWatches && f.$on("$destroy", A
                      .removeWatches)), m) {
                    var E = b[D],
                      O = m[D],
                      N = E.$$bindings.bindToController;
                    O.instance = O(), x.data("$" + E.name + "Controller", O.instance), O.bindingInfo = _e(g, S, O.instance, N, E)
                  }
                  for (C(b, (function(e, t) {
                      var n = e.require;
                      e.bindToController && !G(n) && V(n) && M(m[t].instance, be(t, n, x, m))
                    })), C(m, (function(e) {
                      var t = e.instance;
                      if (z(t.$onChanges)) try {
                        t.$onChanges(e.bindingInfo.initialChanges)
                      } catch (e) {
                        h(e)
                      }
                      if (z(t.$onInit)) try {
                        t.$onInit()
                      } catch (e) {
                        h(e)
                      }
                      z(t.$doCheck) && (g.$watch((function() {
                        t.$doCheck()
                      })), t.$doCheck()), z(t.$onDestroy) && g.$on("$destroy", (function() {
                        t.$onDestroy()
                      }))
                    })), c = 0, d = u.length; c < d; c++) Re(p = u[c], p.isolateScope ? f : t, x, S, p.require && be(p.directiveName, p.require, x, m),
                  v);
                  var P = t;
                  for (y && (y.template || null === y.templateUrl) && (P = f), e && e(P, i.childNodes, void 0, s), c = l.length - 1; c >= 0; c--) Re(p =
                    l[c], p.isolateScope ? f : t, x, S, p.require && be(p.directiveName, p.require, x, m), v);
                  C(m, (function(e) {
                    var t = e.instance;
                    z(t.$postLink) && t.$postLink()
                  }))
                }
              }

              function be(e, t, n, r) {
                var i;
                if (H(t)) {
                  var o = t.match(c),
                    a = t.substring(o[0].length),
                    s = o[1] || o[3],
                    u = "?" === o[2];
                  if ("^^" === s ? n = n.parent() : i = (i = r && r[a]) && i.instance, !i) {
                    var l = "$" + a + "Controller";
                    i = "^^" === s && n[0] && 9 === n[0].nodeType ? null : s ? n.inheritedData(l) : n.data(l)
                  }
                  if (!i && !u) throw gn("ctreq", "Controller '{0}', required by directive '{1}', can't be found!", a, e)
                } else if (G(t)) {
                  i = [];
                  for (var d = 0, p = t.length; d < p; d++) i[d] = be(e, t[d], n, r)
                } else V(t) && (i = {}, C(t, (function(t, o) {
                  i[o] = be(e, t, n, r)
                })));
                return i || null
              }

              function ye(e, t, n) {
                for (var r = 0, i = e.length; r < i; r++) e[r] = P(e[r], {
                  $$isolateScope: t,
                  $$newScope: n
                })
              }

              function xe(e, n, o, a, s, u, l) {
                if (n === s) return null;
                var c = null;
                if (r.hasOwnProperty(n))
                  for (var d, p = t.get(n + i), f = 0, h = p.length; f < h; f++)
                    if (d = p[f], (j(a) || a > d.priority) && -1 !== d.restrict.indexOf(o)) {
                      if (u && (d = P(d, {
                          $$start: u,
                          $$end: l
                        })), !d.$$bindings) {
                        var g = d.$$bindings = m(d, d.name);
                        V(g.isolateScope) && (d.$$isolateBindings = g.isolateScope)
                      }
                      e.push(d), c = d
                    } return c
              }

              function Ce(e) {
                if (r.hasOwnProperty(e))
                  for (var n = t.get(e + i), o = 0, a = n.length; o < a; o++)
                    if (n[o].multiElement) return !0;
                return !1
              }

              function Se(e, t) {
                var n = t.$attr,
                  r = e.$attr;
                C(e, (function(r, i) {
                  "$" !== i.charAt(0) && (t[i] && t[i] !== r && (r.length ? r += ("style" === i ? ";" : " ") + t[i] : r = t[i]), e.$set(i, r, !0, n[
                    i]))
                })), C(t, (function(t, i) {
                  e.hasOwnProperty(i) || "$" === i.charAt(0) || (e[i] = t, "class" !== i && "style" !== i && (r[i] = n[i]))
                }))
              }

              function Ae(e, t, n, r, i, o, s, u) {
                var l, c, d = [],
                  p = t[0],
                  f = e.shift(),
                  m = P(f, {
                    templateUrl: null,
                    transclude: null,
                    replace: null,
                    $$originalDirective: f
                  }),
                  v = z(f.templateUrl) ? f.templateUrl(t, n) : f.templateUrl,
                  $ = f.templateNamespace;
                return t.empty(), g(v).then((function(h) {
                    var g, b, y, w;
                    if (h = te(h), f.replace) {
                      if (y = lt(h) ? [] : Cn(De($, ee(h))), g = y[0], 1 !== y.length || 1 !== g.nodeType) throw gn("tplrt",
                        "Template for directive '{0}' must have exactly one root element. {1}", f.name, v);
                      b = {
                        $attr: {}
                      }, qe(r, t, g);
                      var x = he(g, [], b);
                      V(f.scope) && ye(x, !0), e = x.concat(e), Se(n, b)
                    } else g = p, t.html(h);
                    for (e.unshift(m), l = $e(e, g, n, i, t, f, o, s, u), C(r, (function(e, n) {
                        e === g && (r[n] = t[0])
                      })), c = le(t[0].childNodes, i); d.length;) {
                      var S = d.shift(),
                        A = d.shift(),
                        T = d.shift(),
                        k = d.shift(),
                        D = t[0];
                      if (!S.$$destroyed) {
                        if (A !== p) {
                          var M = A.className;
                          u.hasElementTranscludeDirective && f.replace || (D = ht(g)), qe(T, a(A), D), K(a(D), M)
                        }
                        w = l.transcludeOnThisElement ? de(S, l.transclude, k) : k, l(c, S, D, r, w)
                      }
                    }
                    d = null
                  })).catch((function(e) {
                    W(e) && h(e)
                  })),
                  function(e, t, n, r, i) {
                    var o = i;
                    t.$$destroyed || (d ? d.push(t, n, r, o) : (l.transcludeOnThisElement && (o = de(t, l.transclude, i)), l(c, t, n, r, o)))
                  }
              }

              function Te(e, t) {
                var n = t.priority - e.priority;
                return 0 !== n ? n : e.name !== t.name ? e.name < t.name ? -1 : 1 : e.index - t.index
              }

              function ke(e, t, n, r) {
                function i(e) {
                  return e ? " (module: " + e + ")" : ""
                }
                if (t) throw gn("multidir", "Multiple directives [{0}{1}, {2}{3}] asking for {4} on: {5}", t.name, i(t.$$moduleName), n.name, i(n
                  .$$moduleName), e, we(r))
              }

              function De(t, n) {
                switch (t = p(t || "html")) {
                  case "svg":
                  case "math":
                    var r = e.document.createElement("div");
                    return r.innerHTML = "<" + t + ">" + n + "</" + t + ">", r.childNodes[0].childNodes;
                  default:
                    return n
                }
              }

              function Me(e) {
                return B(D.valueOf(e), "ng-prop-srcset")
              }

              function Ee(e, t, n, r) {
                if (f.test(r)) throw gn("nodomevents", "Property bindings for HTML DOM event properties are disallowed");
                var i = re(e),
                  o = function(e, t) {
                    var n = t.toLowerCase();
                    return S[e + "|" + n] || S["*|" + n]
                  }(i, r),
                  a = I;
                "srcset" !== r || "img" !== i && "source" !== i ? o && (a = D.getTrusted.bind(D, o)) : a = Me, t.push({
                  priority: 100,
                  compile: function(e, t) {
                    var i = A(t[n]),
                      o = A(t[n], (function(e) {
                        return D.valueOf(e)
                      }));
                    return {
                      pre: function(e, t) {
                        function n() {
                          var n = i(e);
                          t[0][r] = a(n)
                        }
                        n(), e.$watch(o, n)
                      }
                    }
                  }
                })
              }

              function Oe(e, t, n) {
                e.push(Ro(A, k, h, t, n, !1))
              }

              function Ne(e, t, r, i, o) {
                var a = re(e),
                  s = function(e, t) {
                    return "srcdoc" === t ? D.HTML : "src" === t || "ngSrc" === t ? -1 === ["img", "video", "audio", "source", "track"].indexOf(e) ? D
                      .RESOURCE_URL : D.MEDIA_URL : "xlinkHref" === t ? "image" === e ? D.MEDIA_URL : "a" === e ? D.URL : D.RESOURCE_URL : "form" === e &&
                      "action" === t || "base" === e && "href" === t || "link" === e && "href" === t ? D.RESOURCE_URL : "a" !== e || "href" !== t &&
                      "ngHref" !== t ? void 0 : D.URL
                  }(a, i),
                  u = !o,
                  c = l[i] || o,
                  d = n(r, u, s, c);
                if (d) {
                  if ("multiple" === i && "select" === a) throw gn("selmulti", "Binding to the 'multiple' attribute is not supported. Element: {0}", we(
                    e));
                  if (f.test(i)) throw gn("nodomevents", "Interpolations for HTML DOM event attributes are disallowed");
                  t.push({
                    priority: 100,
                    compile: function() {
                      return {
                        pre: function(e, t, o) {
                          var a = o.$$observers || (o.$$observers = Ve()),
                            u = o[i];
                          u !== r && (d = u && n(u, !0, s, c), r = u), d && (o[i] = d(e), (a[i] || (a[i] = [])).$$inter = !0, (o.$$observers && o
                            .$$observers[i].$$scope || e).$watch(d, (function(e, t) {
                            "class" === i && e !== t ? o.$updateClass(e, t) : o.$set(i, e)
                          })))
                        }
                      }
                    }
                  })
                }
              }

              function qe(t, n, r) {
                var i, o, s = n[0],
                  u = n.length,
                  l = s.parentNode;
                if (t)
                  for (i = 0, o = t.length; i < o; i++)
                    if (t[i] === s) {
                      t[i++] = r;
                      for (var c = i, d = c + u - 1, p = t.length; c < p; c++, d++) d < p ? t[c] = t[d] : delete t[c];
                      t.length -= u - 1, t.context === s && (t.context = r);
                      break
                    } l && l.replaceChild(r, s);
                var f = e.document.createDocumentFragment();
                for (i = 0; i < u; i++) f.appendChild(n[i]);
                for (a.hasData(s) && (a.data(r, a.data(s)), a(s).off("$destroy")), a.cleanData(f.querySelectorAll("*")), i = 1; i < u; i++) delete n[i];
                n[0] = r, n.length = 1
              }

              function Ie(e, t) {
                return M((function() {
                  return e.apply(null, arguments)
                }), e, t)
              }

              function Re(e, t, n, r, i, o) {
                try {
                  e(t, n, r, i, o)
                } catch (e) {
                  h(e, we(n))
                }
              }

              function je(e, t) {
                if (b) throw gn("missingattr", "Attribute '{0}' of '{1}' is non-optional and must be set!", e, t)
              }

              function _e(e, t, r, i, o) {
                var a, s = [],
                  u = {};

                function l(t, n, i) {
                  z(r.$onChanges) && !se(n, i) && (O || (e.$$postDigest(F), O = []), a || (a = {}, O.push(c)), a[t] && (i = a[t].previousValue), a[t] =
                    new $n(i, n))
                }

                function c() {
                  r.$onChanges(a), a = void 0
                }
                return C(i, (function(i, a) {
                  var c, p, f, h, g, m = i.attrName,
                    v = i.optional;
                  switch (i.mode) {
                    case "@":
                      v || d.call(t, m) || (je(m, o.name), r[a] = t[m] = void 0), g = t.$observe(m, (function(e) {
                        if (H(e) || X(e)) {
                          var t = r[a];
                          l(a, e, t), r[a] = e
                        }
                      })), t.$$observers[m].$$scope = e, H(c = t[m]) ? r[a] = n(c)(e) : X(c) && (r[a] = c), u[a] = new $n(mn, r[a]), s.push(g);
                      break;
                    case "=":
                      if (!d.call(t, m)) {
                        if (v) break;
                        je(m, o.name), t[m] = void 0
                      }
                      if (v && !t[m]) break;
                      p = A(t[m]), h = p.literal ? ue : se, f = p.assign || function() {
                        throw c = r[a] = p(e), gn("nonassign", "Expression '{0}' in attribute '{1}' used with directive '{2}' is non-assignable!",
                          t[m], m, o.name)
                      }, c = r[a] = p(e);
                      var $ = function(t) {
                        return h(t, r[a]) || (h(t, c) ? f(e, t = r[a]) : r[a] = t), c = t
                      };
                      $.$stateful = !0, g = i.collection ? e.$watchCollection(t[m], $) : e.$watch(A(t[m], $), null, p.literal), s.push(g);
                      break;
                    case "<":
                      if (!d.call(t, m)) {
                        if (v) break;
                        je(m, o.name), t[m] = void 0
                      }
                      if (v && !t[m]) break;
                      var b = (p = A(t[m])).literal,
                        y = r[a] = p(e);
                      u[a] = new $n(mn, r[a]), g = e[i.collection ? "$watchCollection" : "$watch"](p, (function(e, t) {
                        if (t === e) {
                          if (t === y || b && ue(t, y)) return;
                          t = y
                        }
                        l(a, e, t), r[a] = e
                      })), s.push(g);
                      break;
                    case "&":
                      if (v || d.call(t, m) || je(m, o.name), (p = t.hasOwnProperty(m) ? A(t[m]) : q) === q && v) break;
                      r[a] = function(t) {
                        return p(e, t)
                      }
                  }
                })), {
                  initialChanges: u,
                  removeWatches: s.length && function() {
                    for (var e = 0, t = s.length; e < t; ++e) s[e]()
                  }
                }
              }
            }
          ]
      }

      function $n(e, t) {
        this.previousValue = e, this.currentValue = t
      }
      vn.$inject = ["$provide", "$$sanitizeUriProvider"], $n.prototype.isFirstChange = function() {
        return this.previousValue === mn
      };
      var bn = /^((?:x|data)[:\-_])/i,
        yn = /[:\-_]+(.)/g;

      function wn(e) {
        return e.replace(bn, "").replace(yn, (function(e, t, n) {
          return n ? t.toUpperCase() : t
        }))
      }

      function xn(e, t) {
        var n = "",
          r = e.split(/\s+/),
          i = t.split(/\s+/);
        e: for (var o = 0; o < r.length; o++) {
          for (var a = r[o], s = 0; s < i.length; s++)
            if (a === i[s]) continue e;
          n += (n.length > 0 ? " " : "") + a
        }
        return n
      }

      function Cn(e) {
        var t = (e = a(e)).length;
        if (t <= 1) return e;
        for (; t--;) {
          var n = e[t];
          (8 === n.nodeType || n.nodeType === He && "" === n.nodeValue.trim()) && g.call(e, t, 1)
        }
        return e
      }
      var Sn = i("$controller"),
        An = /^(\S+)(\s+as\s+([\w$]+))?$/;

      function Tn(e, t) {
        if (t && H(t)) return t;
        if (H(e)) {
          var n = An.exec(e);
          if (n) return n[3]
        }
      }

      function kn() {
        var e = {};
        this.has = function(t) {
          return e.hasOwnProperty(t)
        }, this.register = function(t, n) {
          je(t, "controller"), V(t) ? M(e, t) : e[t] = n
        }, this.$get = ["$injector", function(t) {
          return function(r, i, o, a) {
            var s, u, l, c;
            if (o = !0 === o, a && H(a) && (c = a), H(r)) {
              if (!(u = r.match(An))) throw Sn("ctrlfmt", "Badly formed controller string '{0}'. Must match `__name__ as __id__` or `__name__`.",
              r);
              if (l = u[1], c = c || u[3], !(r = e.hasOwnProperty(l) ? e[l] : function(e, t, n) {
                  if (!t) return e;
                  for (var r, i = t.split("."), o = e, a = i.length, s = 0; s < a; s++) r = i[s], e && (e = (o = e)[r]);
                  return !n && z(e) ? fe(o, e) : e
                }(i.$scope, l, !0))) throw Sn("ctrlreg", "The controller with the name '{0}' is not registered.", l);
              Re(r, l, !0)
            }
            if (o) {
              var d = (G(r) ? r[r.length - 1] : r).prototype;
              return s = Object.create(d || null), c && n(i, c, s, l || r.name), M((function() {
                var e = t.invoke(r, s, i, l);
                return e !== s && (V(e) || z(e)) && (s = e, c && n(i, c, s, l || r.name)), s
              }), {
                instance: s,
                identifier: c
              })
            }
            return s = t.instantiate(r, i, l), c && n(i, c, s, l || r.name), s
          };

          function n(e, t, n, r) {
            if (!e || !V(e.$scope)) throw i("$controller")("noscp",
              "Cannot export controller '{0}' as '{1}'! No $scope object provided via `locals`.", r, t);
            e.$scope[t] = n
          }
        }]
      }

      function Dn() {
        this.$get = ["$window", function(e) {
          return a(e.document)
        }]
      }

      function Mn() {
        this.$get = ["$document", "$rootScope", function(e, t) {
          var n = e[0],
            r = n && n.hidden;

          function i() {
            r = n.hidden
          }
          return e.on("visibilitychange", i), t.$on("$destroy", (function() {
              e.off("visibilitychange", i)
            })),
            function() {
              return r
            }
        }]
      }

      function En() {
        this.$get = ["$log", function(e) {
          return function(t, n) {
            e.error.apply(e, arguments)
          }
        }]
      }
      var On = function() {
          this.$get = ["$document", function(e) {
            return function(t) {
              return t ? !t.nodeType && t instanceof a && (t = t[0]) : t = e[0].body, t.offsetWidth + 1
            }
          }]
        },
        Nn = "application/json",
        Pn = {
          "Content-Type": "application/json;charset=utf-8"
        },
        qn = /^\[|^\{(?!\{)/,
        In = {
          "[": /]$/,
          "{": /}$/
        },
        Ln = /^\)]\}',?\n/,
        Rn = i("$http");

      function jn(e) {
        return V(e) ? B(e) ? e.toISOString() : ge(e) : e
      }

      function _n() {
        this.$get = function() {
          return function(e) {
            if (!e) return "";
            var t = [];
            return S(e, (function(e, n) {
              null === e || j(e) || z(e) || (G(e) ? C(e, (function(e) {
                t.push(Ae(n) + "=" + Ae(jn(e)))
              })) : t.push(Ae(n) + "=" + Ae(jn(e))))
            })), t.join("&")
          }
        }
      }

      function Vn() {
        this.$get = function() {
          return function(e) {
            if (!e) return "";
            var t = [];
            return function e(n, r, i) {
              G(n) ? C(n, (function(t, n) {
                e(t, r + "[" + (V(t) ? n : "") + "]")
              })) : V(n) && !B(n) ? S(n, (function(t, n) {
                e(t, r + (i ? "" : "[") + n + (i ? "" : "]"))
              })) : (z(n) && (n = n()), t.push(Ae(r) + "=" + (null == n ? "" : Ae(jn(n)))))
            }(e, "", !0), t.join("&")
          }
        }
      }

      function Un(e, t) {
        if (H(e)) {
          var n = e.replace(Ln, "").trim();
          if (n) {
            var r = t("Content-Type"),
              i = r && 0 === r.indexOf(Nn);
            if (i || (a = (o = n).match(qn)) && In[a[0]].test(o)) try {
              e = me(n)
            } catch (t) {
              if (!i) return e;
              throw Rn("baddata", 'Data must be a valid JSON object. Received: "{0}". Parse error: "{1}"', e, t)
            }
          }
        }
        var o, a;
        return e
      }

      function Hn(e) {
        var t, n = Ve();

        function r(e, t) {
          e && (n[e] = n[e] ? n[e] + ", " + t : t)
        }
        return H(e) ? C(e.split("\n"), (function(e) {
          t = e.indexOf(":"), r(p(ee(e.substr(0, t))), ee(e.substr(t + 1)))
        })) : V(e) && C(e, (function(e, t) {
          r(p(t), ee(e))
        })), n
      }

      function Fn(e) {
        var t;
        return function(n) {
          if (t || (t = Hn(e)), n) {
            var r = t[p(n)];
            return void 0 === r && (r = null), r
          }
          return t
        }
      }

      function Bn(e, t, n, r) {
        return z(r) ? r(e, t, n) : (C(r, (function(r) {
          e = r(e, t, n)
        })), e)
      }

      function Gn(e) {
        return 200 <= e && e < 300
      }

      function Wn() {
        var e = this.defaults = {
            transformResponse: [Un],
            transformRequest: [function(e) {
              return !V(e) || (t = e, "[object File]" === v.call(t)) || function(e) {
                return "[object Blob]" === v.call(e)
              }(e) || function(e) {
                return "[object FormData]" === v.call(e)
              }(e) ? e : ge(e);
              var t
            }],
            headers: {
              common: {
                Accept: "application/json, text/plain, */*"
              },
              post: Fe(Pn),
              put: Fe(Pn),
              patch: Fe(Pn)
            },
            xsrfCookieName: "XSRF-TOKEN",
            xsrfHeaderName: "X-XSRF-TOKEN",
            paramSerializer: "$httpParamSerializer",
            jsonpCallbackParam: "callback"
          },
          t = !1;
        this.useApplyAsync = function(e) {
          return _(e) ? (t = !!e, this) : t
        };
        var n = this.interceptors = [],
          r = this.xsrfTrustedOrigins = [];
        Object.defineProperty(this, "xsrfWhitelistedOrigins", {
          get: function() {
            return this.xsrfTrustedOrigins
          },
          set: function(e) {
            this.xsrfTrustedOrigins = e
          }
        }), this.$get = ["$browser", "$httpBackend", "$$cookieReader", "$cacheFactory", "$rootScope", "$q", "$injector", "$sce", function(o, a, s, u, l,
          c, d, h) {
          var g = u("$http");
          e.paramSerializer = H(e.paramSerializer) ? d.get(e.paramSerializer) : e.paramSerializer;
          var m = [];
          C(n, (function(e) {
            m.unshift(H(e) ? d.get(e) : d.invoke(e))
          }));
          var v, $ = (v = [li].concat(r.map(di)), function(e) {
            var t = di(e);
            return v.some(pi.bind(null, t))
          });

          function b(n) {
            if (!V(n)) throw i("$http")("badreq", "Http request configuration must be an object.  Received: {0}", n);
            if (!H(h.valueOf(n.url))) throw i("$http")("badreq",
              "Http request configuration url must be a string or a $sce trusted object.  Received: {0}", n.url);
            var r = M({
              method: "get",
              transformRequest: e.transformRequest,
              transformResponse: e.transformResponse,
              paramSerializer: e.paramSerializer,
              jsonpCallbackParam: e.jsonpCallbackParam
            }, n);
            r.headers = function(t) {
                var n, r, i, o = e.headers,
                  a = M({}, t.headers);
                o = M({}, o.common, o[p(t.method)]);
                e: for (n in o) {
                  for (i in r = p(n), a)
                    if (p(i) === r) continue e;
                  a[n] = o[n]
                }
                return function(e, t) {
                  var n, r = {};
                  return C(e, (function(e, i) {
                    z(e) ? null != (n = e(t)) && (r[i] = n) : r[i] = e
                  })), r
                }(a, Fe(t))
              }(n), r.method = f(r.method), r.paramSerializer = H(r.paramSerializer) ? d.get(r.paramSerializer) : r.paramSerializer, o
              .$$incOutstandingRequestCount("$http");
            var u = [],
              v = [],
              y = c.resolve(r);
            return C(m, (function(e) {
              (e.request || e.requestError) && u.unshift(e.request, e.requestError), (e.response || e.responseError) && v.push(e.response, e
                .responseError)
            })), y = w(y, u), y = (y = w(y = y.then((function(n) {
              var r = n.headers,
                i = Bn(n.data, Fn(r), void 0, n.transformRequest);
              j(i) && C(r, (function(e, t) {
                "content-type" === p(t) && delete r[t]
              }));
              j(n.withCredentials) && !j(e.withCredentials) && (n.withCredentials = e.withCredentials);
              return function(n, r) {
                var i, o, u = c.defer(),
                  d = u.promise,
                  f = n.headers,
                  m = "jsonp" === p(n.method),
                  v = n.url;
                m ? v = h.getTrustedResourceUrl(v) : H(v) || (v = h.valueOf(v));
                v = function(e, t) {
                  t.length > 0 && (e += (-1 === e.indexOf("?") ? "?" : "&") + t);
                  return e
                }(v, n.paramSerializer(n.params)), m && (v = function(e, t) {
                  var n = e.split("?");
                  if (n.length > 2) throw Rn("badjsonp", 'Illegal use more than one "?", in url, "{1}"', e);
                  return C(Ce(n[1]), (function(n, r) {
                    if ("JSON_CALLBACK" === n) throw Rn("badjsonp", 'Illegal use of JSON_CALLBACK in url, "{0}"', e);
                    if (r === t) throw Rn("badjsonp", 'Illegal use of callback param, "{0}", in url, "{1}"', t, e)
                  })), e += (-1 === e.indexOf("?") ? "?" : "&") + t + "=JSON_CALLBACK"
                }(v, n.jsonpCallbackParam));
                b.pendingRequests.push(n), d.then(T, T), !n.cache && !e.cache || !1 === n.cache || "GET" !== n.method && "JSONP" !== n
                  .method || (i = V(n.cache) ? n.cache : V(e.cache) ? e.cache : g);
                i && (_(o = i.get(v)) ? J(o) ? o.then(A, A) : G(o) ? S(o[1], o[0], Fe(o[2]), o[3], o[4]) : S(o, 200, {}, "OK",
                  "complete") : i.put(v, d));
                if (j(o)) {
                  var y = $(n.url) ? s()[n.xsrfCookieName || e.xsrfCookieName] : void 0;
                  y && (f[n.xsrfHeaderName || e.xsrfHeaderName] = y), a(n.method, v, r, x, f, n.timeout, n.withCredentials, n
                    .responseType, w(n.eventHandlers), w(n.uploadEventHandlers))
                }
                return d;

                function w(e) {
                  if (e) {
                    var n = {};
                    return C(e, (function(e, r) {
                      n[r] = function(n) {
                        function r() {
                          e(n)
                        }
                        t ? l.$applyAsync(r) : l.$$phase ? r() : l.$apply(r)
                      }
                    })), n
                  }
                }

                function x(e, n, r, o, a) {
                  function s() {
                    S(n, e, r, o, a)
                  }
                  i && (Gn(e) ? i.put(v, [e, n, Hn(r), o, a]) : i.remove(v)), t ? l.$applyAsync(s) : (s(), l.$$phase || l.$apply())
                }

                function S(e, t, r, i, o) {
                  (Gn(t = t >= -1 ? t : 0) ? u.resolve : u.reject)({
                    data: e,
                    status: t,
                    headers: Fn(r),
                    config: n,
                    statusText: i,
                    xhrStatus: o
                  })
                }

                function A(e) {
                  S(e.data, e.status, Fe(e.headers()), e.statusText, e.xhrStatus)
                }

                function T() {
                  var e = b.pendingRequests.indexOf(n); - 1 !== e && b.pendingRequests.splice(e, 1)
                }
              }(n, i).then(x, x)
            })), v)).finally((function() {
              o.$$completeOutstandingRequest(q, "$http")
            }));

            function w(e, t) {
              for (var n = 0, r = t.length; n < r;) {
                var i = t[n++],
                  o = t[n++];
                e = e.then(i, o)
              }
              return t.length = 0, e
            }

            function x(e) {
              var t = M({}, e);
              return t.data = Bn(e.data, e.headers, e.status, r.transformResponse), Gn(e.status) ? t : c.reject(t)
            }
          }
          return b.pendingRequests = [],
            function(e) {
              C(arguments, (function(e) {
                b[e] = function(t, n) {
                  return b(M({}, n || {}, {
                    method: e,
                    url: t
                  }))
                }
              }))
            }("get", "delete", "head", "jsonp"),
            function(e) {
              C(arguments, (function(e) {
                b[e] = function(t, n, r) {
                  return b(M({}, r || {}, {
                    method: e,
                    url: t,
                    data: n
                  }))
                }
              }))
            }("post", "put", "patch"), b.defaults = e, b
        }]
      }

      function zn() {
        this.$get = function() {
          return function() {
            return new e.XMLHttpRequest
          }
        }
      }

      function Yn() {
        this.$get = ["$browser", "$jsonpCallbacks", "$document", "$xhrFactory", function(e, t, n, r) {
          return function(e, t, n, r, i) {
            return function(i, a, s, u, l, c, d, f, h, g) {
              if (a = a || e.url(), "jsonp" === p(i)) var m = r.createCallback(a),
                v = o(a, m, (function(e, t) {
                  var n = 200 === e && r.getResponse(m);
                  T(u, e, n, "", t, "complete"), r.removeCallback(m)
                }));
              else {
                var $ = t(i, a),
                  b = !1;
                $.open(i, a, !0), C(l, (function(e, t) {
                  _(e) && $.setRequestHeader(t, e)
                })), $.onload = function() {
                  var e = $.statusText || "",
                    t = "response" in $ ? $.response : $.responseText,
                    n = 1223 === $.status ? 204 : $.status;
                  0 === n && (n = t ? 200 : "file" === di(a).protocol ? 404 : 0), T(u, n, t, $.getAllResponseHeaders(), e, "complete")
                };
                var y = function() {
                    T(u, -1, null, null, "", "error")
                  },
                  w = function() {
                    T(u, -1, null, null, "", b ? "timeout" : "abort")
                  },
                  x = function() {
                    T(u, -1, null, null, "", "timeout")
                  };
                if ($.onerror = y, $.ontimeout = x, $.onabort = w, C(h, (function(e, t) {
                    $.addEventListener(t, e)
                  })), C(g, (function(e, t) {
                    $.upload.addEventListener(t, e)
                  })), d && ($.withCredentials = !0), f) try {
                  $.responseType = f
                } catch (e) {
                  if ("json" !== f) throw e
                }
                $.send(j(s) ? null : s)
              }
              if (c > 0) var S = n((function() {
                A("timeout")
              }), c);
              else J(c) && c.then((function() {
                A(_(c.$$timeoutId) ? "timeout" : "abort")
              }));

              function A(e) {
                b = "timeout" === e, v && v(), $ && $.abort()
              }

              function T(e, t, r, i, o, a) {
                _(S) && n.cancel(S), v = $ = null, e(t, r, i, o, a)
              }
            };

            function o(e, t, n) {
              e = e.replace("JSON_CALLBACK", t);
              var o = i.createElement("script"),
                a = null;
              return o.type = "text/javascript", o.src = e, o.async = !0, a = function(e) {
                o.removeEventListener("load", a), o.removeEventListener("error", a), i.body.removeChild(o), o = null;
                var s = -1,
                  u = "unknown";
                e && ("load" !== e.type || r.wasCalled(t) || (e = {
                  type: "error"
                }), u = e.type, s = "error" === e.type ? 404 : 200), n && n(s, u)
              }, o.addEventListener("load", a), o.addEventListener("error", a), i.body.appendChild(o), a
            }
          }(e, r, e.defer, t, n[0])
        }]
      }
      var Kn = y.$interpolateMinErr = i("$interpolate");

      function Qn() {
        var e = "{{",
          t = "}}";
        this.startSymbol = function(t) {
          return t ? (e = t, this) : e
        }, this.endSymbol = function(e) {
          return e ? (t = e, this) : t
        }, this.$get = ["$parse", "$exceptionHandler", "$sce", function(n, r, i) {
          var o = e.length,
            a = t.length,
            s = new RegExp(e.replace(/./g, l), "g"),
            u = new RegExp(t.replace(/./g, l), "g");

          function l(e) {
            return "\\\\\\" + e
          }

          function c(n) {
            return n.replace(s, e).replace(u, t)
          }

          function d(e, t, n, r) {
            var i = e.$watch((function(e) {
              return i(), r(e)
            }), t, n);
            return i
          }

          function p(s, u, l, p) {
            var f = l === i.URL || l === i.MEDIA_URL;
            if (!s.length || -1 === s.indexOf(e)) {
              if (u) return;
              var h = c(s);
              f && (h = i.getTrusted(l, h));
              var g = L(h);
              return g.exp = s, g.expressions = [], g.$$watchDelegate = d, g
            }
            p = !!p;
            for (var m, v, $, b, y, w = 0, x = [], C = s.length, S = [], A = []; w < C;) {
              if (-1 === (m = s.indexOf(e, w)) || -1 === (v = s.indexOf(t, m + o))) {
                w !== C && S.push(c(s.substring(w)));
                break
              }
              w !== m && S.push(c(s.substring(w, m))), b = s.substring(m + o, v), x.push(b), w = v + a, A.push(S.length), S.push("")
            }
            y = 1 === S.length && 1 === A.length;
            var T = f && y ? void 0 : function(e) {
              try {
                return e = l && !f ? i.getTrusted(l, e) : i.valueOf(e), p && !_(e) ? e : Ue(e)
              } catch (e) {
                r(Kn.interr(s, e))
              }
            };
            if ($ = x.map((function(e) {
                return n(e, T)
              })), !u || x.length) {
              var k = function(e) {
                for (var t = 0, n = x.length; t < n; t++) {
                  if (p && j(e[t])) return;
                  S[A[t]] = e[t]
                }
                return f ? i.getTrusted(l, y ? S[0] : S.join("")) : (l && S.length > 1 && Kn.throwNoconcat(s), S.join(""))
              };
              return M((function(e) {
                var t = 0,
                  n = x.length,
                  i = new Array(n);
                try {
                  for (; t < n; t++) i[t] = $[t](e);
                  return k(i)
                } catch (e) {
                  r(Kn.interr(s, e))
                }
              }), {
                exp: s,
                expressions: x,
                $$watchDelegate: function(e, t) {
                  var n;
                  return e.$watchGroup($, (function(r, i) {
                    var o = k(r);
                    t.call(this, o, r !== i ? n : o, e), n = o
                  }))
                }
              })
            }
          }
          return p.startSymbol = function() {
            return e
          }, p.endSymbol = function() {
            return t
          }, p
        }]
      }
      Kn.throwNoconcat = function(e) {
        throw Kn("noconcat",
          "Error while interpolating: {0}\nStrict Contextual Escaping disallows interpolations that concatenate multiple expressions when a trusted value is required.  See http://docs.angularjs.org/api/ng.$sce",
          e)
      }, Kn.interr = function(e, t) {
        return Kn("interr", "Can't interpolate: {0}\n{1}", e, t.toString())
      };
      var Xn = i("$interval");

      function Jn() {
        this.$get = ["$$intervalFactory", "$window", function(e, t) {
          var n = {},
            r = function(e) {
              t.clearInterval(e), delete n[e]
            },
            i = e((function(e, r, i) {
              var o = t.setInterval(e, r);
              return n[o] = i, o
            }), r);
          return i.cancel = function(e) {
            if (!e) return !1;
            if (!e.hasOwnProperty("$$intervalId")) throw Xn("badprom",
              "`$interval.cancel()` called with a promise that was not generated by `$interval()`.");
            if (!n.hasOwnProperty(e.$$intervalId)) return !1;
            var t = e.$$intervalId,
              i = n[t];
            return Hr(i.promise), i.reject("canceled"), r(t), !0
          }, i
        }]
      }

      function Zn() {
        this.$get = ["$browser", "$q", "$$q", "$rootScope", function(e, t, n, r) {
          return function(i, o) {
            return function(a, s, u, l) {
              var c = arguments.length > 4,
                d = c ? pe(arguments, 4) : [],
                p = 0,
                f = _(l) && !l,
                h = (f ? n : t).defer(),
                g = h.promise;

              function m() {
                c ? a.apply(null, d) : a(p)
              }

              function v() {
                f ? e.defer(m) : r.$evalAsync(m), h.notify(p++), u > 0 && p >= u && (h.resolve(p), o(g.$$intervalId)), f || r.$apply()
              }
              return u = _(u) ? u : 0, g.$$intervalId = i(v, s, h, f), g
            }
          }
        }]
      }
      var er = function() {
          this.$get = function() {
            var e = y.callbacks,
              t = {};
            return {
              createCallback: function(n) {
                var r = "_" + (e.$$counter++).toString(36),
                  i = "angular.callbacks." + r,
                  o = function(e) {
                    var t = function(e) {
                      t.data = e, t.called = !0
                    };
                    return t.id = e, t
                  }(r);
                return t[i] = e[r] = o, i
              },
              wasCalled: function(e) {
                return t[e].called
              },
              getResponse: function(e) {
                return t[e].data
              },
              removeCallback: function(n) {
                var r = t[n];
                delete e[r.id], delete t[n]
              }
            }
          }
        },
        tr = /^([^?#]*)(\?([^#]*))?(#(.*))?$/,
        nr = {
          http: 80,
          https: 443,
          ftp: 21
        },
        rr = i("$location");

      function ir(e, t, n) {
        var r, i = (r = [], C(t, (function(e, t) {
            G(e) ? C(e, (function(e) {
              r.push(Ae(t, !0) + (!0 === e ? "" : "=" + Ae(e, !0)))
            })) : r.push(Ae(t, !0) + (!0 === e ? "" : "=" + Ae(e, !0)))
          })), r.length ? r.join("&") : ""),
          o = n ? "#" + Se(n) : "";
        return function(e) {
          for (var t = e.split("/"), n = t.length; n--;) t[n] = Se(t[n].replace(/%2F/g, "/"));
          return t.join("/")
        }(e) + (i ? "?" + i : "") + o
      }

      function or(e, t) {
        var n = di(e);
        t.$$protocol = n.protocol, t.$$host = n.hostname, t.$$port = O(n.port) || nr[n.protocol] || null
      }
      var ar = /^\s*[\\/]{2,}/;

      function sr(e, t, n) {
        if (ar.test(e)) throw rr("badpath", 'Invalid url "{0}".', e);
        var r = "/" !== e.charAt(0);
        r && (e = "/" + e);
        var i = di(e),
          o = r && "/" === i.pathname.charAt(0) ? i.pathname.substring(1) : i.pathname;
        t.$$path = function(e, t) {
          for (var n = e.split("/"), r = n.length; r--;) n[r] = decodeURIComponent(n[r]), t && (n[r] = n[r].replace(/\//g, "%2F"));
          return n.join("/")
        }(o, n), t.$$search = Ce(i.search), t.$$hash = decodeURIComponent(i.hash), t.$$path && "/" !== t.$$path.charAt(0) && (t.$$path = "/" + t.$$path)
      }

      function ur(e, t) {
        return e.slice(0, t.length) === t
      }

      function lr(e, t) {
        if (ur(t, e)) return t.substr(e.length)
      }

      function cr(e) {
        var t = e.indexOf("#");
        return -1 === t ? e : e.substr(0, t)
      }

      function dr(e, t, n) {
        this.$$html5 = !0, n = n || "", or(e, this), this.$$parse = function(e) {
          var n = lr(t, e);
          if (!H(n)) throw rr("ipthprfx", 'Invalid url "{0}", missing path prefix "{1}".', e, t);
          sr(n, this, !0), this.$$path || (this.$$path = "/"), this.$$compose()
        }, this.$$normalizeUrl = function(e) {
          return t + e.substr(1)
        }, this.$$parseLinkUrl = function(r, i) {
          return i && "#" === i[0] ? (this.hash(i.slice(1)), !0) : (_(o = lr(e, r)) ? (a = o, s = n && _(o = lr(n, o)) ? t + (lr("/", o) || o) : e +
            a) : _(o = lr(t, r)) ? s = t + o : t === r + "/" && (s = t), s && this.$$parse(s), !!s);
          var o, a, s
        }
      }

      function pr(e, t, n) {
        or(e, this), this.$$parse = function(r) {
          var i, o = lr(e, r) || lr(t, r);
          j(o) || "#" !== o.charAt(0) ? this.$$html5 ? i = o : (i = "", j(o) && (e = r, this.replace())) : j(i = lr(n, o)) && (i = o), sr(i, this, !1),
            this.$$path = function(e, t, n) {
              var r, i = /^\/[A-Z]:(\/.*)/;
              ur(t, n) && (t = t.replace(n, ""));
              if (i.exec(t)) return e;
              return (r = i.exec(e)) ? r[1] : e
            }(this.$$path, i, e), this.$$compose()
        }, this.$$normalizeUrl = function(t) {
          return e + (t ? n + t : "")
        }, this.$$parseLinkUrl = function(t, n) {
          return cr(e) === cr(t) && (this.$$parse(t), !0)
        }
      }

      function fr(e, t, n) {
        this.$$html5 = !0, pr.apply(this, arguments), this.$$parseLinkUrl = function(r, i) {
          return i && "#" === i[0] ? (this.hash(i.slice(1)), !0) : (e === cr(r) ? o = r : (a = lr(t, r)) ? o = e + n + a : t === r + "/" && (o = t),
            o && this.$$parse(o), !!o);
          var o, a
        }, this.$$normalizeUrl = function(t) {
          return e + n + t
        }
      }
      var hr = {
        $$absUrl: "",
        $$html5: !1,
        $$replace: !1,
        $$compose: function() {
          this.$$url = ir(this.$$path, this.$$search, this.$$hash), this.$$absUrl = this.$$normalizeUrl(this.$$url), this.$$urlUpdatedByLocation = !0
        },
        absUrl: gr("$$absUrl"),
        url: function(e) {
          if (j(e)) return this.$$url;
          var t = tr.exec(e);
          return (t[1] || "" === e) && this.path(decodeURIComponent(t[1])), (t[2] || t[1] || "" === e) && this.search(t[3] || ""), this.hash(t[5] ||
            ""), this
        },
        protocol: gr("$$protocol"),
        host: gr("$$host"),
        port: gr("$$port"),
        path: mr("$$path", (function(e) {
          return "/" === (e = null !== e ? e.toString() : "").charAt(0) ? e : "/" + e
        })),
        search: function(e, t) {
          switch (arguments.length) {
            case 0:
              return this.$$search;
            case 1:
              if (H(e) || F(e)) e = e.toString(), this.$$search = Ce(e);
              else {
                if (!V(e)) throw rr("isrcharg", "The first argument of the `$location#search()` call must be a string or an object.");
                C(e = ae(e, {}), (function(t, n) {
                  null == t && delete e[n]
                })), this.$$search = e
              }
              break;
            default:
              j(t) || null === t ? delete this.$$search[e] : this.$$search[e] = t
          }
          return this.$$compose(), this
        },
        hash: mr("$$hash", (function(e) {
          return null !== e ? e.toString() : ""
        })),
        replace: function() {
          return this.$$replace = !0, this
        }
      };

      function gr(e) {
        return function() {
          return this[e]
        }
      }

      function mr(e, t) {
        return function(n) {
          return j(n) ? this[e] : (this[e] = t(n), this.$$compose(), this)
        }
      }

      function vr() {
        var e = "!",
          t = {
            enabled: !1,
            requireBase: !0,
            rewriteLinks: !0
          };
        this.hashPrefix = function(t) {
          return _(t) ? (e = t, this) : e
        }, this.html5Mode = function(e) {
          return X(e) ? (t.enabled = e, this) : V(e) ? (X(e.enabled) && (t.enabled = e.enabled), X(e.requireBase) && (t.requireBase = e.requireBase), (
            X(e.rewriteLinks) || H(e.rewriteLinks)) && (t.rewriteLinks = e.rewriteLinks), this) : t
        }, this.$get = ["$rootScope", "$browser", "$sniffer", "$rootElement", "$window", function(n, r, i, o, s) {
          var u, l, c, d, p = r.baseHref(),
            f = r.url();
          if (t.enabled) {
            if (!p && t.requireBase) throw rr("nobase", "$location in HTML5 mode requires a <base> tag to be present!");
            c = (d = f).substring(0, d.indexOf("/", d.indexOf("//") + 2)) + (p || "/"), l = i.history ? dr : fr
          } else c = cr(f), l = pr;
          var h = function(e) {
            return e.substr(0, cr(e).lastIndexOf("/") + 1)
          }(c);
          (u = new l(c, h, "#" + e)).$$parseLinkUrl(f, f), u.$$state = r.state();
          var g = /^\s*(javascript|mailto):/i;

          function m(e, t, n) {
            var i = u.url(),
              o = u.$$state;
            try {
              r.url(e, t, n), u.$$state = r.state()
            } catch (e) {
              throw u.url(i), u.$$state = o, e
            }
          }
          o.on("click", (function(e) {
            var i = t.rewriteLinks;
            if (i && !e.ctrlKey && !e.metaKey && !e.shiftKey && 2 !== e.which && 2 !== e.button) {
              for (var s = a(e.target);
                "a" !== re(s[0]);)
                if (s[0] === o[0] || !(s = s.parent())[0]) return;
              if (!H(i) || !j(s.attr(i))) {
                var l = s.prop("href"),
                  c = s.attr("href") || s.attr("xlink:href");
                V(l) && "[object SVGAnimatedString]" === l.toString() && (l = di(l.animVal).href), g.test(l) || !l || s.attr("target") || e
                  .isDefaultPrevented() || u.$$parseLinkUrl(l, c) && (e.preventDefault(), u.absUrl() !== r.url() && n.$apply())
              }
            }
          })), u.absUrl() !== f && r.url(u.absUrl(), !0);
          var v = !0;
          return r.onUrlChange((function(e, t) {
            ur(e, h) ? (n.$evalAsync((function() {
              var r, i = u.absUrl(),
                o = u.$$state;
              u.$$parse(e), u.$$state = t, r = n.$broadcast("$locationChangeStart", e, i, t, o).defaultPrevented, u.absUrl() === e && (
                r ? (u.$$parse(i), u.$$state = o, m(i, !1, o)) : (v = !1, $(i, o)))
            })), n.$$phase || n.$digest()) : s.location.href = e
          })), n.$watch((function() {
            if (v || u.$$urlUpdatedByLocation) {
              u.$$urlUpdatedByLocation = !1;
              var e = r.url(),
                t = u.absUrl(),
                o = r.state(),
                a = u.$$replace,
                s = !((l = e) === (c = t) || di(l).href === di(c).href) || u.$$html5 && i.history && o !== u.$$state;
              (v || s) && (v = !1, n.$evalAsync((function() {
                var t = u.absUrl(),
                  r = n.$broadcast("$locationChangeStart", t, e, u.$$state, o).defaultPrevented;
                u.absUrl() === t && (r ? (u.$$parse(e), u.$$state = o) : (s && m(t, a, o === u.$$state ? null : u.$$state), $(e, o)))
              })))
            }
            var l, c;
            u.$$replace = !1
          })), u;

          function $(e, t) {
            n.$broadcast("$locationChangeSuccess", u.absUrl(), e, u.$$state, t)
          }
        }]
      }

      function $r() {
        var e = !0,
          t = this;
        this.debugEnabled = function(t) {
          return _(t) ? (e = t, this) : e
        }, this.$get = ["$window", function(n) {
          var r, i = o || /\bEdge\//.test(n.navigator && n.navigator.userAgent);
          return {
            log: s("log"),
            info: s("info"),
            warn: s("warn"),
            error: s("error"),
            debug: (r = s("debug"), function() {
              e && r.apply(t, arguments)
            })
          };

          function a(e) {
            return W(e) && (e.stack && i ? e = e.message && -1 === e.stack.indexOf(e.message) ? "Error: " + e.message + "\n" + e.stack : e.stack : e
              .sourceURL && (e = e.message + "\n" + e.sourceURL + ":" + e.line)), e
          }

          function s(e) {
            var t = n.console || {},
              r = t[e] || t.log || q;
            return function() {
              var e = [];
              return C(arguments, (function(t) {
                e.push(a(t))
              })), Function.prototype.apply.call(r, t, e)
            }
          }
        }]
      }
      C([fr, pr, dr], (function(e) {
        e.prototype = Object.create(hr), e.prototype.state = function(t) {
          if (!arguments.length) return this.$$state;
          if (e !== dr || !this.$$html5) throw rr("nostate",
            "History API state support is available only in HTML5 mode and only in browsers supporting HTML5 History API");
          return this.$$state = j(t) ? null : t, this.$$urlUpdatedByLocation = !0, this
        }
      }));
      var br = i("$parse"),
        yr = {}.constructor.prototype.valueOf;

      function wr(e) {
        return e + ""
      }
      var xr = Ve();
      C("+ - * / % === !== == != < > <= >= && || ! = |".split(" "), (function(e) {
        xr[e] = !0
      }));
      var Cr = {
          n: "\n",
          f: "\f",
          r: "\r",
          t: "\t",
          v: "\v",
          "'": "'",
          '"': '"'
        },
        Sr = function(e) {
          this.options = e
        };
      Sr.prototype = {
        constructor: Sr,
        lex: function(e) {
          for (this.text = e, this.index = 0, this.tokens = []; this.index < this.text.length;) {
            var t = this.text.charAt(this.index);
            if ('"' === t || "'" === t) this.readString(t);
            else if (this.isNumber(t) || "." === t && this.isNumber(this.peek())) this.readNumber();
            else if (this.isIdentifierStart(this.peekMultichar())) this.readIdent();
            else if (this.is(t, "(){}[].,;:?")) this.tokens.push({
              index: this.index,
              text: t
            }), this.index++;
            else if (this.isWhitespace(t)) this.index++;
            else {
              var n = t + this.peek(),
                r = n + this.peek(2),
                i = xr[t],
                o = xr[n],
                a = xr[r];
              if (i || o || a) {
                var s = a ? r : o ? n : t;
                this.tokens.push({
                  index: this.index,
                  text: s,
                  operator: !0
                }), this.index += s.length
              } else this.throwError("Unexpected next character ", this.index, this.index + 1)
            }
          }
          return this.tokens
        },
        is: function(e, t) {
          return -1 !== t.indexOf(e)
        },
        peek: function(e) {
          var t = e || 1;
          return this.index + t < this.text.length && this.text.charAt(this.index + t)
        },
        isNumber: function(e) {
          return "0" <= e && e <= "9" && "string" == typeof e
        },
        isWhitespace: function(e) {
          return " " === e || "\r" === e || "\t" === e || "\n" === e || "\v" === e || " " === e
        },
        isIdentifierStart: function(e) {
          return this.options.isIdentifierStart ? this.options.isIdentifierStart(e, this.codePointAt(e)) : this.isValidIdentifierStart(e)
        },
        isValidIdentifierStart: function(e) {
          return "a" <= e && e <= "z" || "A" <= e && e <= "Z" || "_" === e || "$" === e
        },
        isIdentifierContinue: function(e) {
          return this.options.isIdentifierContinue ? this.options.isIdentifierContinue(e, this.codePointAt(e)) : this.isValidIdentifierContinue(e)
        },
        isValidIdentifierContinue: function(e, t) {
          return this.isValidIdentifierStart(e, t) || this.isNumber(e)
        },
        codePointAt: function(e) {
          return 1 === e.length ? e.charCodeAt(0) : (e.charCodeAt(0) << 10) + e.charCodeAt(1) - 56613888
        },
        peekMultichar: function() {
          var e = this.text.charAt(this.index),
            t = this.peek();
          if (!t) return e;
          var n = e.charCodeAt(0),
            r = t.charCodeAt(0);
          return n >= 55296 && n <= 56319 && r >= 56320 && r <= 57343 ? e + t : e
        },
        isExpOperator: function(e) {
          return "-" === e || "+" === e || this.isNumber(e)
        },
        throwError: function(e, t, n) {
          n = n || this.index;
          var r = _(t) ? "s " + t + "-" + this.index + " [" + this.text.substring(t, n) + "]" : " " + n;
          throw br("lexerr", "Lexer Error: {0} at column{1} in expression [{2}].", e, r, this.text)
        },
        readNumber: function() {
          for (var e = "", t = this.index; this.index < this.text.length;) {
            var n = p(this.text.charAt(this.index));
            if ("." === n || this.isNumber(n)) e += n;
            else {
              var r = this.peek();
              if ("e" === n && this.isExpOperator(r)) e += n;
              else if (this.isExpOperator(n) && r && this.isNumber(r) && "e" === e.charAt(e.length - 1)) e += n;
              else {
                if (!this.isExpOperator(n) || r && this.isNumber(r) || "e" !== e.charAt(e.length - 1)) break;
                this.throwError("Invalid exponent")
              }
            }
            this.index++
          }
          this.tokens.push({
            index: t,
            text: e,
            constant: !0,
            value: Number(e)
          })
        },
        readIdent: function() {
          var e = this.index;
          for (this.index += this.peekMultichar().length; this.index < this.text.length;) {
            var t = this.peekMultichar();
            if (!this.isIdentifierContinue(t)) break;
            this.index += t.length
          }
          this.tokens.push({
            index: e,
            text: this.text.slice(e, this.index),
            identifier: !0
          })
        },
        readString: function(e) {
          var t = this.index;
          this.index++;
          for (var n = "", r = e, i = !1; this.index < this.text.length;) {
            var o = this.text.charAt(this.index);
            if (r += o, i) {
              if ("u" === o) {
                var a = this.text.substring(this.index + 1, this.index + 5);
                a.match(/[\da-f]{4}/i) || this.throwError("Invalid unicode escape [\\u" + a + "]"), this.index += 4, n += String.fromCharCode(
                  parseInt(a, 16))
              } else {
                n += Cr[o] || o
              }
              i = !1
            } else if ("\\" === o) i = !0;
            else {
              if (o === e) return this.index++, void this.tokens.push({
                index: t,
                text: r,
                constant: !0,
                value: n
              });
              n += o
            }
            this.index++
          }
          this.throwError("Unterminated quote", t)
        }
      };
      var Ar = function(e, t) {
        this.lexer = e, this.options = t
      };

      function Tr(e, t) {
        return void 0 !== e ? e : t
      }

      function kr(e, t) {
        return void 0 === e ? t : void 0 === t ? e : e + t
      }
      Ar.Program = "Program", Ar.ExpressionStatement = "ExpressionStatement", Ar.AssignmentExpression = "AssignmentExpression", Ar.ConditionalExpression =
        "ConditionalExpression", Ar.LogicalExpression = "LogicalExpression", Ar.BinaryExpression = "BinaryExpression", Ar.UnaryExpression =
        "UnaryExpression", Ar.CallExpression = "CallExpression", Ar.MemberExpression = "MemberExpression", Ar.Identifier = "Identifier", Ar.Literal =
        "Literal", Ar.ArrayExpression = "ArrayExpression", Ar.Property = "Property", Ar.ObjectExpression = "ObjectExpression", Ar.ThisExpression =
        "ThisExpression", Ar.LocalsExpression = "LocalsExpression", Ar.NGValueParameter = "NGValueParameter", Ar.prototype = {
          ast: function(e) {
            this.text = e, this.tokens = this.lexer.lex(e);
            var t = this.program();
            return 0 !== this.tokens.length && this.throwError("is an unexpected token", this.tokens[0]), t
          },
          program: function() {
            for (var e = [];;)
              if (this.tokens.length > 0 && !this.peek("}", ")", ";", "]") && e.push(this.expressionStatement()), !this.expect(";")) return {
                type: Ar.Program,
                body: e
              }
          },
          expressionStatement: function() {
            return {
              type: Ar.ExpressionStatement,
              expression: this.filterChain()
            }
          },
          filterChain: function() {
            for (var e = this.expression(); this.expect("|");) e = this.filter(e);
            return e
          },
          expression: function() {
            return this.assignment()
          },
          assignment: function() {
            var e = this.ternary();
            if (this.expect("=")) {
              if (!Er(e)) throw br("lval", "Trying to assign a value to a non l-value");
              e = {
                type: Ar.AssignmentExpression,
                left: e,
                right: this.assignment(),
                operator: "="
              }
            }
            return e
          },
          ternary: function() {
            var e, t, n = this.logicalOR();
            return this.expect("?") && (e = this.expression(), this.consume(":")) ? (t = this.expression(), {
              type: Ar.ConditionalExpression,
              test: n,
              alternate: e,
              consequent: t
            }) : n
          },
          logicalOR: function() {
            for (var e = this.logicalAND(); this.expect("||");) e = {
              type: Ar.LogicalExpression,
              operator: "||",
              left: e,
              right: this.logicalAND()
            };
            return e
          },
          logicalAND: function() {
            for (var e = this.equality(); this.expect("&&");) e = {
              type: Ar.LogicalExpression,
              operator: "&&",
              left: e,
              right: this.equality()
            };
            return e
          },
          equality: function() {
            for (var e, t = this.relational(); e = this.expect("==", "!=", "===", "!==");) t = {
              type: Ar.BinaryExpression,
              operator: e.text,
              left: t,
              right: this.relational()
            };
            return t
          },
          relational: function() {
            for (var e, t = this.additive(); e = this.expect("<", ">", "<=", ">=");) t = {
              type: Ar.BinaryExpression,
              operator: e.text,
              left: t,
              right: this.additive()
            };
            return t
          },
          additive: function() {
            for (var e, t = this.multiplicative(); e = this.expect("+", "-");) t = {
              type: Ar.BinaryExpression,
              operator: e.text,
              left: t,
              right: this.multiplicative()
            };
            return t
          },
          multiplicative: function() {
            for (var e, t = this.unary(); e = this.expect("*", "/", "%");) t = {
              type: Ar.BinaryExpression,
              operator: e.text,
              left: t,
              right: this.unary()
            };
            return t
          },
          unary: function() {
            var e;
            return (e = this.expect("+", "-", "!")) ? {
              type: Ar.UnaryExpression,
              operator: e.text,
              prefix: !0,
              argument: this.unary()
            } : this.primary()
          },
          primary: function() {
            var e, t;
            for (this.expect("(") ? (e = this.filterChain(), this.consume(")")) : this.expect("[") ? e = this.arrayDeclaration() : this.expect("{") ?
              e = this.object() : this.selfReferential.hasOwnProperty(this.peek().text) ? e = ae(this.selfReferential[this.consume().text]) : this
              .options.literals.hasOwnProperty(this.peek().text) ? e = {
                type: Ar.Literal,
                value: this.options.literals[this.consume().text]
              } : this.peek().identifier ? e = this.identifier() : this.peek().constant ? e = this.constant() : this.throwError(
                "not a primary expression", this.peek()); t = this.expect("(", "[", ".");) "(" === t.text ? (e = {
              type: Ar.CallExpression,
              callee: e,
              arguments: this.parseArguments()
            }, this.consume(")")) : "[" === t.text ? (e = {
              type: Ar.MemberExpression,
              object: e,
              property: this.expression(),
              computed: !0
            }, this.consume("]")) : "." === t.text ? e = {
              type: Ar.MemberExpression,
              object: e,
              property: this.identifier(),
              computed: !1
            } : this.throwError("IMPOSSIBLE");
            return e
          },
          filter: function(e) {
            for (var t = [e], n = {
                type: Ar.CallExpression,
                callee: this.identifier(),
                arguments: t,
                filter: !0
              }; this.expect(":");) t.push(this.expression());
            return n
          },
          parseArguments: function() {
            var e = [];
            if (")" !== this.peekToken().text)
              do {
                e.push(this.filterChain())
              } while (this.expect(","));
            return e
          },
          identifier: function() {
            var e = this.consume();
            return e.identifier || this.throwError("is not a valid identifier", e), {
              type: Ar.Identifier,
              name: e.text
            }
          },
          constant: function() {
            return {
              type: Ar.Literal,
              value: this.consume().value
            }
          },
          arrayDeclaration: function() {
            var e = [];
            if ("]" !== this.peekToken().text)
              do {
                if (this.peek("]")) break;
                e.push(this.expression())
              } while (this.expect(","));
            return this.consume("]"), {
              type: Ar.ArrayExpression,
              elements: e
            }
          },
          object: function() {
            var e, t = [];
            if ("}" !== this.peekToken().text)
              do {
                if (this.peek("}")) break;
                e = {
                    type: Ar.Property,
                    kind: "init"
                  }, this.peek().constant ? (e.key = this.constant(), e.computed = !1, this.consume(":"), e.value = this.expression()) : this.peek()
                  .identifier ? (e.key = this.identifier(), e.computed = !1, this.peek(":") ? (this.consume(":"), e.value = this.expression()) : e
                    .value = e.key) : this.peek("[") ? (this.consume("["), e.key = this.expression(), this.consume("]"), e.computed = !0, this.consume(
                    ":"), e.value = this.expression()) : this.throwError("invalid key", this.peek()), t.push(e)
              } while (this.expect(","));
            return this.consume("}"), {
              type: Ar.ObjectExpression,
              properties: t
            }
          },
          throwError: function(e, t) {
            throw br("syntax", "Syntax Error: Token '{0}' {1} at column {2} of the expression [{3}] starting at [{4}].", t.text, e, t.index + 1, this
              .text, this.text.substring(t.index))
          },
          consume: function(e) {
            if (0 === this.tokens.length) throw br("ueoe", "Unexpected end of expression: {0}", this.text);
            var t = this.expect(e);
            return t || this.throwError("is unexpected, expecting [" + e + "]", this.peek()), t
          },
          peekToken: function() {
            if (0 === this.tokens.length) throw br("ueoe", "Unexpected end of expression: {0}", this.text);
            return this.tokens[0]
          },
          peek: function(e, t, n, r) {
            return this.peekAhead(0, e, t, n, r)
          },
          peekAhead: function(e, t, n, r, i) {
            if (this.tokens.length > e) {
              var o = this.tokens[e],
                a = o.text;
              if (a === t || a === n || a === r || a === i || !t && !n && !r && !i) return o
            }
            return !1
          },
          expect: function(e, t, n, r) {
            var i = this.peek(e, t, n, r);
            return !!i && (this.tokens.shift(), i)
          },
          selfReferential: {
            this: {
              type: Ar.ThisExpression
            },
            $locals: {
              type: Ar.LocalsExpression
            }
          }
        };

      function Dr(e, t, n) {
        var r, i, o, a = e.isPure = function(e, t) {
          switch (e.type) {
            case Ar.MemberExpression:
              if (e.computed) return !1;
              break;
            case Ar.UnaryExpression:
              return 1;
            case Ar.BinaryExpression:
              return "+" !== e.operator && 1;
            case Ar.CallExpression:
              return !1
          }
          return void 0 === t ? 2 : t
        }(e, n);
        switch (e.type) {
          case Ar.Program:
            r = !0, C(e.body, (function(e) {
              Dr(e.expression, t, a), r = r && e.expression.constant
            })), e.constant = r;
            break;
          case Ar.Literal:
            e.constant = !0, e.toWatch = [];
            break;
          case Ar.UnaryExpression:
            Dr(e.argument, t, a), e.constant = e.argument.constant, e.toWatch = e.argument.toWatch;
            break;
          case Ar.BinaryExpression:
            Dr(e.left, t, a), Dr(e.right, t, a), e.constant = e.left.constant && e.right.constant, e.toWatch = e.left.toWatch.concat(e.right.toWatch);
            break;
          case Ar.LogicalExpression:
            Dr(e.left, t, a), Dr(e.right, t, a), e.constant = e.left.constant && e.right.constant, e.toWatch = e.constant ? [] : [e];
            break;
          case Ar.ConditionalExpression:
            Dr(e.test, t, a), Dr(e.alternate, t, a), Dr(e.consequent, t, a), e.constant = e.test.constant && e.alternate.constant && e.consequent
              .constant, e.toWatch = e.constant ? [] : [e];
            break;
          case Ar.Identifier:
            e.constant = !1, e.toWatch = [e];
            break;
          case Ar.MemberExpression:
            Dr(e.object, t, a), e.computed && Dr(e.property, t, a), e.constant = e.object.constant && (!e.computed || e.property.constant), e.toWatch = e
              .constant ? [] : [e];
            break;
          case Ar.CallExpression:
            o = !!e.filter && function(e, t) {
              return !e(t).$stateful
            }(t, e.callee.name), r = o, i = [], C(e.arguments, (function(e) {
              Dr(e, t, a), r = r && e.constant, i.push.apply(i, e.toWatch)
            })), e.constant = r, e.toWatch = o ? i : [e];
            break;
          case Ar.AssignmentExpression:
            Dr(e.left, t, a), Dr(e.right, t, a), e.constant = e.left.constant && e.right.constant, e.toWatch = [e];
            break;
          case Ar.ArrayExpression:
            r = !0, i = [], C(e.elements, (function(e) {
              Dr(e, t, a), r = r && e.constant, i.push.apply(i, e.toWatch)
            })), e.constant = r, e.toWatch = i;
            break;
          case Ar.ObjectExpression:
            r = !0, i = [], C(e.properties, (function(e) {
              Dr(e.value, t, a), r = r && e.value.constant, i.push.apply(i, e.value.toWatch), e.computed && (Dr(e.key, t, !1), r = r && e.key
                .constant, i.push.apply(i, e.key.toWatch))
            })), e.constant = r, e.toWatch = i;
            break;
          case Ar.ThisExpression:
          case Ar.LocalsExpression:
            e.constant = !1, e.toWatch = []
        }
      }

      function Mr(e) {
        if (1 === e.length) {
          var t = e[0].expression,
            n = t.toWatch;
          return 1 !== n.length || n[0] !== t ? n : void 0
        }
      }

      function Er(e) {
        return e.type === Ar.Identifier || e.type === Ar.MemberExpression
      }

      function Or(e) {
        if (1 === e.body.length && Er(e.body[0].expression)) return {
          type: Ar.AssignmentExpression,
          left: e.body[0].expression,
          right: {
            type: Ar.NGValueParameter
          },
          operator: "="
        }
      }

      function Nr(e) {
        this.$filter = e
      }

      function Pr(e) {
        this.$filter = e
      }

      function qr(e, t, n) {
        this.ast = new Ar(e, n), this.astCompiler = n.csp ? new Pr(t) : new Nr(t)
      }

      function Ir(e) {
        return z(e.valueOf) ? e.valueOf() : yr.call(e)
      }

      function Lr() {
        var e, t, n = Ve(),
          r = {
            true: !0,
            false: !1,
            null: null,
            undefined: void 0
          };
        this.addLiteral = function(e, t) {
          r[e] = t
        }, this.setIdentifierFns = function(n, r) {
          return e = n, t = r, this
        }, this.$get = ["$filter", function(i) {
          var o = {
            csp: le().noUnsafeEval,
            literals: ae(r),
            isIdentifierStart: z(e) && e,
            isIdentifierContinue: z(t) && t
          };
          return a.$$getAst = function(e) {
            return new qr(new Sr(o), i, o).getAst(e).ast
          }, a;

          function a(e, t) {
            var r, a;
            switch (typeof e) {
              case "string":
                if (e = e.trim(), !(r = n[a = e])) r = new qr(new Sr(o), i, o).parse(e), n[a] = p(r);
                return f(r, t);
              case "function":
                return f(e, t);
              default:
                return f(q, t)
            }
          }

          function s(e, t, n) {
            return null == e || null == t ? e === t : !("object" == typeof e && "object" == typeof(e = Ir(e)) && !n) && (e === t || e != e && t != t)
          }

          function u(e, t, n, r, i) {
            var o, a = r.inputs;
            if (1 === a.length) {
              var u = s;
              return a = a[0], e.$watch((function(e) {
                var t = a(e);
                return s(t, u, a.isPure) || (o = r(e, void 0, void 0, [t]), u = t && Ir(t)), o
              }), t, n, i)
            }
            for (var l = [], c = [], d = 0, p = a.length; d < p; d++) l[d] = s, c[d] = null;
            return e.$watch((function(e) {
              for (var t = !1, n = 0, i = a.length; n < i; n++) {
                var u = a[n](e);
                (t || (t = !s(u, l[n], a[n].isPure))) && (c[n] = u, l[n] = u && Ir(u))
              }
              return t && (o = r(e, void 0, void 0, c)), o
            }), t, n, i)
          }

          function l(e, t, n, r, i) {
            var o, a, s = r.literal ? c : _,
              u = r.$$intercepted || r,
              l = r.$$interceptor || I,
              d = r.inputs && !u.inputs;
            return h.literal = r.literal, h.constant = r.constant, h.inputs = r.inputs, p(h), o = e.$watch(h, t, n, i);

            function f() {
              s(a) && o()
            }

            function h(e, t, n, r) {
              return a = d && r ? r[0] : u(e, t, n, r), s(a) && e.$$postDigest(f), l(a)
            }
          }

          function c(e) {
            var t = !0;
            return C(e, (function(e) {
              _(e) || (t = !1)
            })), t
          }

          function d(e, t, n, r) {
            var i = e.$watch((function(e) {
              return i(), r(e)
            }), t, n);
            return i
          }

          function p(e) {
            return e.constant ? e.$$watchDelegate = d : e.oneTime ? e.$$watchDelegate = l : e.inputs && (e.$$watchDelegate = u), e
          }

          function f(e, t) {
            if (!t) return e;
            e.$$interceptor && (t = function(e, t) {
              function n(n) {
                return t(e(n))
              }
              return n.$stateful = e.$stateful || t.$stateful, n.$$pure = e.$$pure && t.$$pure, n
            }(e.$$interceptor, t), e = e.$$intercepted);
            var n = !1,
              r = function(r, i, o, a) {
                var s = n && a ? a[0] : e(r, i, o, a);
                return t(s)
              };
            return r.$$intercepted = e, r.$$interceptor = t, r.literal = e.literal, r.oneTime = e.oneTime, r.constant = e.constant, t.$stateful || (
              n = !e.inputs, r.inputs = e.inputs ? e.inputs : [e], t.$$pure || (r.inputs = r.inputs.map((function(e) {
                return 2 === e.isPure ? function(t) {
                  return e(t)
                } : e
              })))), p(r)
          }
        }]
      }

      function Rr() {
        var e = !0;
        this.$get = ["$rootScope", "$exceptionHandler", function(t, n) {
          return _r((function(e) {
            t.$evalAsync(e)
          }), n, e)
        }], this.errorOnUnhandledRejections = function(t) {
          return _(t) ? (e = t, this) : e
        }
      }

      function jr() {
        var e = !0;
        this.$get = ["$browser", "$exceptionHandler", function(t, n) {
          return _r((function(e) {
            t.defer(e)
          }), n, e)
        }], this.errorOnUnhandledRejections = function(t) {
          return _(t) ? (e = t, this) : e
        }
      }

      function _r(e, t, n) {
        var r = i("$q", TypeError),
          o = 0,
          a = [];

        function s() {
          return new u
        }

        function u() {
          var e = this.promise = new l;
          this.resolve = function(t) {
            p(e, t)
          }, this.reject = function(t) {
            h(e, t)
          }, this.notify = function(t) {
            m(e, t)
          }
        }

        function l() {
          this.$$state = {
            status: 0
          }
        }

        function c() {
          for (; !o && a.length;) {
            var e = a.shift();
            if (!Vr(e)) {
              Ur(e);
              var n = "Possibly unhandled rejection: " + Be(e.value);
              W(e.value) ? t(e.value, n) : t(n)
            }
          }
        }

        function d(r) {
          !n || r.pending || 2 !== r.status || Vr(r) || (0 === o && 0 === a.length && e(c), a.push(r)), !r.processScheduled && r.pending && (r
            .processScheduled = !0, ++o, e((function() {
              ! function(r) {
                var i, a, s;
                s = r.pending, r.processScheduled = !1, r.pending = void 0;
                try {
                  for (var u = 0, l = s.length; u < l; ++u) {
                    Ur(r), a = s[u][0], i = s[u][r.status];
                    try {
                      z(i) ? p(a, i(r.value)) : 1 === r.status ? p(a, r.value) : h(a, r.value)
                    } catch (e) {
                      h(a, e), e && !0 === e.$$passToExceptionHandler && t(e)
                    }
                  }
                } finally {
                  --o, n && 0 === o && e(c)
                }
              }(r)
            })))
        }

        function p(e, t) {
          e.$$state.status || (t === e ? g(e, r("qcycle", "Expected promise to be resolved with value other than itself '{0}'", t)) : f(e, t))
        }

        function f(e, t) {
          var n, r = !1;
          try {
            (V(t) || z(t)) && (n = t.then), z(n) ? (e.$$state.status = -1, n.call(t, (function(t) {
              if (r) return;
              r = !0, f(e, t)
            }), i, (function(t) {
              m(e, t)
            }))) : (e.$$state.value = t, e.$$state.status = 1, d(e.$$state))
          } catch (e) {
            i(e)
          }

          function i(t) {
            r || (r = !0, g(e, t))
          }
        }

        function h(e, t) {
          e.$$state.status || g(e, t)
        }

        function g(e, t) {
          e.$$state.value = t, e.$$state.status = 2, d(e.$$state)
        }

        function m(n, r) {
          var i = n.$$state.pending;
          n.$$state.status <= 0 && i && i.length && e((function() {
            for (var e, n, o = 0, a = i.length; o < a; o++) {
              n = i[o][0], e = i[o][3];
              try {
                m(n, z(e) ? e(r) : r)
              } catch (e) {
                t(e)
              }
            }
          }))
        }

        function v(e) {
          var t = new l;
          return h(t, e), t
        }

        function $(e, t, n) {
          var r = null;
          try {
            z(n) && (r = n())
          } catch (e) {
            return v(e)
          }
          return J(r) ? r.then((function() {
            return t(e)
          }), v) : t(e)
        }

        function b(e, t, n, r) {
          var i = new l;
          return p(i, e), i.then(t, n, r)
        }
        M(l.prototype, {
          then: function(e, t, n) {
            if (j(e) && j(t) && j(n)) return this;
            var r = new l;
            return this.$$state.pending = this.$$state.pending || [], this.$$state.pending.push([r, e, t, n]), this.$$state.status > 0 && d(this
              .$$state), r
          },
          catch: function(e) {
            return this.then(null, e)
          },
          finally: function(e, t) {
            return this.then((function(t) {
              return $(t, y, e)
            }), (function(t) {
              return $(t, v, e)
            }), t)
          }
        });
        var y = b;

        function w(e) {
          if (!z(e)) throw r("norslvr", "Expected resolverFn, got '{0}'", e);
          var t = new l;
          return e((function(e) {
            p(t, e)
          }), (function(e) {
            h(t, e)
          })), t
        }
        return w.prototype = l.prototype, w.defer = s, w.reject = v, w.when = b, w.resolve = y, w.all = function(e) {
          var t = new l,
            n = 0,
            r = G(e) ? [] : {};
          return C(e, (function(e, i) {
            n++, b(e).then((function(e) {
              r[i] = e, --n || p(t, r)
            }), (function(e) {
              h(t, e)
            }))
          })), 0 === n && p(t, r), t
        }, w.race = function(e) {
          var t = s();
          return C(e, (function(e) {
            b(e).then(t.resolve, t.reject)
          })), t.promise
        }, w
      }

      function Vr(e) {
        return !!e.pur
      }

      function Ur(e) {
        e.pur = !0
      }

      function Hr(e) {
        e.$$state && Ur(e.$$state)
      }

      function Fr() {
        this.$get = ["$window", "$timeout", function(e, t) {
          var n = e.requestAnimationFrame || e.webkitRequestAnimationFrame,
            r = e.cancelAnimationFrame || e.webkitCancelAnimationFrame || e.webkitCancelRequestAnimationFrame,
            i = !!n,
            o = i ? function(e) {
              var t = n(e);
              return function() {
                r(t)
              }
            } : function(e) {
              var n = t(e, 16.66, !1);
              return function() {
                t.cancel(n)
              }
            };
          return o.supported = i, o
        }]
      }

      function Br() {
        var e = 10,
          t = i("$rootScope"),
          n = null,
          r = null;
        this.digestTtl = function(t) {
          return arguments.length && (e = t), e
        }, this.$get = ["$exceptionHandler", "$parse", "$browser", function(i, a, s) {
          function u(e) {
            e.currentScope.$$destroyed = !0
          }

          function l(e) {
            9 === o && (e.$$childHead && l(e.$$childHead), e.$$nextSibling && l(e.$$nextSibling)), e.$parent = e.$$nextSibling = e.$$prevSibling = e
              .$$childHead = e.$$childTail = e.$root = e.$$watchers = null
          }

          function c() {
            this.$id = T(), this.$$phase = this.$parent = this.$$watchers = this.$$nextSibling = this.$$prevSibling = this.$$childHead = this
              .$$childTail = null, this.$root = this, this.$$destroyed = !1, this.$$suspended = !1, this.$$listeners = {}, this.$$listenerCount = {},
              this.$$watchersCount = 0, this.$$isolateBindings = null
          }
          c.prototype = {
            constructor: c,
            $new: function(e, t) {
              var n;
              return t = t || this, e ? (n = new c).$root = this.$root : (this.$$ChildScope || (this.$$ChildScope = function(e) {
                function t() {
                  this.$$watchers = this.$$nextSibling = this.$$childHead = this.$$childTail = null, this.$$listeners = {}, this
                    .$$listenerCount = {}, this.$$watchersCount = 0, this.$id = T(), this.$$ChildScope = null, this.$$suspended = !1
                }
                return t.prototype = e, t
              }(this)), n = new this.$$ChildScope), n.$parent = t, n.$$prevSibling = t.$$childTail, t.$$childHead ? (t.$$childTail
                .$$nextSibling = n, t.$$childTail = n) : t.$$childHead = t.$$childTail = n, (e || t !== this) && n.$on("$destroy", u), n
            },
            $watch: function(e, t, r, i) {
              var o = a(e),
                s = z(t) ? t : q;
              if (o.$$watchDelegate) return o.$$watchDelegate(this, s, r, o, e);
              var u = this,
                l = u.$$watchers,
                c = {
                  fn: s,
                  last: w,
                  get: o,
                  exp: i || e,
                  eq: !!r
                };
              return n = null, l || ((l = u.$$watchers = []).$$digestWatchIndex = -1), l.unshift(c), l.$$digestWatchIndex++, b(this, 1),
                function() {
                  var e = oe(l, c);
                  e >= 0 && (b(u, -1), e < l.$$digestWatchIndex && l.$$digestWatchIndex--), n = null
                }
            },
            $watchGroup: function(e, t) {
              var n = new Array(e.length),
                r = new Array(e.length),
                i = [],
                o = this,
                a = !1,
                s = !0;
              if (!e.length) {
                var u = !0;
                return o.$evalAsync((function() {
                    u && t(r, r, o)
                  })),
                  function() {
                    u = !1
                  }
              }
              if (1 === e.length) return this.$watch(e[0], (function(e, i, o) {
                r[0] = e, n[0] = i, t(r, e === i ? r : n, o)
              }));

              function l() {
                a = !1;
                try {
                  s ? (s = !1, t(r, r, o)) : t(r, n, o)
                } finally {
                  for (var i = 0; i < e.length; i++) n[i] = r[i]
                }
              }
              return C(e, (function(e, t) {
                  var n = o.$watch(e, (function(e) {
                    r[t] = e, a || (a = !0, o.$evalAsync(l))
                  }));
                  i.push(n)
                })),
                function() {
                  for (; i.length;) i.shift()()
                }
            },
            $watchCollection: function(e, t) {
              g.$$pure = a(e).literal, g.$stateful = !g.$$pure;
              var n, r, i, o = this,
                s = t.length > 1,
                u = 0,
                l = a(e, g),
                c = [],
                p = {},
                f = !0,
                h = 0;

              function g(e) {
                var t, i, o, a;
                if (!j(n = e)) {
                  if (V(n))
                    if (x(n)) {
                      r !== c && (h = (r = c).length = 0, u++), t = n.length, h !== t && (u++, r.length = h = t);
                      for (var s = 0; s < t; s++) a = r[s], o = n[s], a != a && o != o || a === o || (u++, r[s] = o)
                    } else {
                      for (i in r !== p && (r = p = {}, h = 0, u++), t = 0, n) d.call(n, i) && (t++, o = n[i], a = r[i], i in r ? a != a && o !=
                        o || a === o || (u++, r[i] = o) : (h++, r[i] = o, u++));
                      if (h > t)
                        for (i in u++, r) d.call(n, i) || (h--, delete r[i])
                    }
                  else r !== n && (r = n, u++);
                  return u
                }
              }
              return this.$watch(l, (function() {
                if (f ? (f = !1, t(n, n, o)) : t(n, i, o), s)
                  if (V(n))
                    if (x(n)) {
                      i = new Array(n.length);
                      for (var e = 0; e < n.length; e++) i[e] = n[e]
                    } else
                      for (var r in i = {}, n) d.call(n, r) && (i[r] = n[r]);
                else i = n
              }))
            },
            $digest: function() {
              var o, a, u, l, c, d, g, b, y, x = e,
                C = f.length ? p : this,
                A = [];
              v("$digest"), s.$$checkUrlChange(), this === p && null !== r && (s.defer.cancel(r), S()), n = null;
              do {
                c = !1, g = C;
                for (var T = 0; T < f.length; T++) {
                  try {
                    (0, (y = f[T]).fn)(y.scope, y.locals)
                  } catch (e) {
                    i(e)
                  }
                  n = null
                }
                f.length = 0;
                e: do {
                  if (l = !g.$$suspended && g.$$watchers)
                    for (l.$$digestWatchIndex = l.length; l.$$digestWatchIndex--;) try {
                      if (o = l[l.$$digestWatchIndex])
                        if ((a = (0, o.get)(g)) === (u = o.last) || (o.eq ? ue(a, u) : N(a) && N(u))) {
                          if (o === n) {
                            c = !1;
                            break e
                          }
                        } else c = !0, n = o, o.last = o.eq ? ae(a, null) : a, (0, o.fn)(a, u === w ? a : u, g), x < 5 && (A[b = 4 - x] || (A[
                          b] = []), A[b].push({
                          msg: z(o.exp) ? "fn: " + (o.exp.name || o.exp.toString()) : o.exp,
                          newVal: a,
                          oldVal: u
                        }))
                    } catch (e) {
                      i(e)
                    }
                  if (!(d = !g.$$suspended && g.$$watchersCount && g.$$childHead || g !== C && g.$$nextSibling))
                    for (; g !== C && !(d = g.$$nextSibling);) g = g.$parent
                } while (g = d);
                if ((c || f.length) && !x--) throw $(), t("infdig",
                  "{0} $digest() iterations reached. Aborting!\nWatchers fired in the last 5 iterations: {1}", e, A)
              } while (c || f.length);
              for ($(); m < h.length;) try {
                h[m++]()
              } catch (e) {
                i(e)
              }
              h.length = m = 0, s.$$checkUrlChange()
            },
            $suspend: function() {
              this.$$suspended = !0
            },
            $isSuspended: function() {
              return this.$$suspended
            },
            $resume: function() {
              this.$$suspended = !1
            },
            $destroy: function() {
              if (!this.$$destroyed) {
                var e = this.$parent;
                for (var t in this.$broadcast("$destroy"), this.$$destroyed = !0, this === p && s.$$applicationDestroyed(), b(this, -this
                    .$$watchersCount), this.$$listenerCount) y(this, this.$$listenerCount[t], t);
                e && e.$$childHead === this && (e.$$childHead = this.$$nextSibling), e && e.$$childTail === this && (e.$$childTail = this
                    .$$prevSibling), this.$$prevSibling && (this.$$prevSibling.$$nextSibling = this.$$nextSibling), this.$$nextSibling && (this
                    .$$nextSibling.$$prevSibling = this.$$prevSibling), this.$destroy = this.$digest = this.$apply = this.$evalAsync = this
                  .$applyAsync = q, this.$on = this.$watch = this.$watchGroup = function() {
                    return q
                  }, this.$$listeners = {}, this.$$nextSibling = null, l(this)
              }
            },
            $eval: function(e, t) {
              return a(e)(this, t)
            },
            $evalAsync: function(e, t) {
              p.$$phase || f.length || s.defer((function() {
                f.length && p.$digest()
              }), null, "$evalAsync"), f.push({
                scope: this,
                fn: a(e),
                locals: t
              })
            },
            $$postDigest: function(e) {
              h.push(e)
            },
            $apply: function(e) {
              try {
                v("$apply");
                try {
                  return this.$eval(e)
                } finally {
                  $()
                }
              } catch (e) {
                i(e)
              } finally {
                try {
                  p.$digest()
                } catch (e) {
                  throw i(e), e
                }
              }
            },
            $applyAsync: function(e) {
              var t = this;
              e && g.push((function() {
                t.$eval(e)
              })), e = a(e), null === r && (r = s.defer((function() {
                p.$apply(S)
              }), null, "$applyAsync"))
            },
            $on: function(e, t) {
              var n = this.$$listeners[e];
              n || (this.$$listeners[e] = n = []), n.push(t);
              var r = this;
              do {
                r.$$listenerCount[e] || (r.$$listenerCount[e] = 0), r.$$listenerCount[e]++
              } while (r = r.$parent);
              var i = this;
              return function() {
                var r = n.indexOf(t); - 1 !== r && (delete n[r], y(i, 1, e))
              }
            },
            $emit: function(e, t) {
              var n, r, o, a = [],
                s = this,
                u = !1,
                l = {
                  name: e,
                  targetScope: s,
                  stopPropagation: function() {
                    u = !0
                  },
                  preventDefault: function() {
                    l.defaultPrevented = !0
                  },
                  defaultPrevented: !1
                },
                c = de([l], arguments, 1);
              do {
                for (n = s.$$listeners[e] || a, l.currentScope = s, r = 0, o = n.length; r < o; r++)
                  if (n[r]) try {
                    n[r].apply(null, c)
                  } catch (e) {
                    i(e)
                  } else n.splice(r, 1), r--, o--;
                if (u) break;
                s = s.$parent
              } while (s);
              return l.currentScope = null, l
            },
            $broadcast: function(e, t) {
              var n = this,
                r = n,
                o = n,
                a = {
                  name: e,
                  targetScope: n,
                  preventDefault: function() {
                    a.defaultPrevented = !0
                  },
                  defaultPrevented: !1
                };
              if (!n.$$listenerCount[e]) return a;
              for (var s, u, l, c = de([a], arguments, 1); r = o;) {
                for (a.currentScope = r, u = 0, l = (s = r.$$listeners[e] || []).length; u < l; u++)
                  if (s[u]) try {
                    s[u].apply(null, c)
                  } catch (e) {
                    i(e)
                  } else s.splice(u, 1), u--, l--;
                if (!(o = r.$$listenerCount[e] && r.$$childHead || r !== n && r.$$nextSibling))
                  for (; r !== n && !(o = r.$$nextSibling);) r = r.$parent
              }
              return a.currentScope = null, a
            }
          };
          var p = new c,
            f = p.$$asyncQueue = [],
            h = p.$$postDigestQueue = [],
            g = p.$$applyAsyncQueue = [],
            m = 0;
          return p;

          function v(e) {
            if (p.$$phase) throw t("inprog", "{0} already in progress", p.$$phase);
            p.$$phase = e
          }

          function $() {
            p.$$phase = null
          }

          function b(e, t) {
            do {
              e.$$watchersCount += t
            } while (e = e.$parent)
          }

          function y(e, t, n) {
            do {
              e.$$listenerCount[n] -= t, 0 === e.$$listenerCount[n] && delete e.$$listenerCount[n]
            } while (e = e.$parent)
          }

          function w() {}

          function S() {
            for (; g.length;) try {
              g.shift()()
            } catch (e) {
              i(e)
            }
            r = null
          }
        }]
      }

      function Gr() {
        var e = /^\s*(https?|s?ftp|mailto|tel|file):/,
          t = /^\s*((https?|ftp|file|blob):|data:image\/)/;
        this.aHrefSanitizationTrustedUrlList = function(t) {
          return _(t) ? (e = t, this) : e
        }, this.imgSrcSanitizationTrustedUrlList = function(e) {
          return _(e) ? (t = e, this) : t
        }, this.$get = function() {
          return function(n, r) {
            var i = r ? t : e,
              o = di(n && n.trim()).href;
            return "" === o || o.match(i) ? n : "unsafe:" + o
          }
        }
      }
      Nr.prototype = {
        compile: function(e) {
          var t = this;
          this.state = {
            nextId: 0,
            filters: {},
            fn: {
              vars: [],
              body: [],
              own: {}
            },
            assign: {
              vars: [],
              body: [],
              own: {}
            },
            inputs: []
          }, Dr(e, t.$filter);
          var n, r = "";
          if (this.stage = "assign", n = Or(e)) {
            this.state.computing = "assign";
            var i = this.nextId();
            this.recurse(n, i), this.return_(i), r = "fn.assign=" + this.generateFunction("assign", "s,v,l")
          }
          var o = Mr(e.body);
          t.stage = "inputs", C(o, (function(e, n) {
            var r = "fn" + n;
            t.state[r] = {
              vars: [],
              body: [],
              own: {}
            }, t.state.computing = r;
            var i = t.nextId();
            t.recurse(e, i), t.return_(i), t.state.inputs.push({
              name: r,
              isPure: e.isPure
            }), e.watchId = n
          })), this.state.computing = "fn", this.stage = "main", this.recurse(e);
          var a = '"' + this.USE + " " + this.STRICT + '";\n' + this.filterPrefix() + "var fn=" + this.generateFunction("fn", "s,l,a,i") + r + this
            .watchFns() + "return fn;",
            s = new Function("$filter", "getStringValue", "ifDefined", "plus", a)(this.$filter, wr, Tr, kr);
          return this.state = this.stage = void 0, s
        },
        USE: "use",
        STRICT: "strict",
        watchFns: function() {
          var e = [],
            t = this.state.inputs,
            n = this;
          return C(t, (function(t) {
            e.push("var " + t.name + "=" + n.generateFunction(t.name, "s")), t.isPure && e.push(t.name, ".isPure=" + JSON.stringify(t.isPure) +
              ";")
          })), t.length && e.push("fn.inputs=[" + t.map((function(e) {
            return e.name
          })).join(",") + "];"), e.join("")
        },
        generateFunction: function(e, t) {
          return "function(" + t + "){" + this.varsPrefix(e) + this.body(e) + "};"
        },
        filterPrefix: function() {
          var e = [],
            t = this;
          return C(this.state.filters, (function(n, r) {
            e.push(n + "=$filter(" + t.escape(r) + ")")
          })), e.length ? "var " + e.join(",") + ";" : ""
        },
        varsPrefix: function(e) {
          return this.state[e].vars.length ? "var " + this.state[e].vars.join(",") + ";" : ""
        },
        body: function(e) {
          return this.state[e].body.join("")
        },
        recurse: function(e, t, n, r, i, o) {
          var a, s, u, l, c, d = this;
          if (r = r || q, !o && _(e.watchId)) return t = t || this.nextId(), void this.if_("i", this.lazyAssign(t, this.computedMember("i", e
            .watchId)), this.lazyRecurse(e, t, n, r, i, !0));
          switch (e.type) {
            case Ar.Program:
              C(e.body, (function(t, n) {
                d.recurse(t.expression, void 0, void 0, (function(e) {
                  s = e
                })), n !== e.body.length - 1 ? d.current().body.push(s, ";") : d.return_(s)
              }));
              break;
            case Ar.Literal:
              l = this.escape(e.value), this.assign(t, l), r(t || l);
              break;
            case Ar.UnaryExpression:
              this.recurse(e.argument, void 0, void 0, (function(e) {
                s = e
              })), l = e.operator + "(" + this.ifDefined(s, 0) + ")", this.assign(t, l), r(l);
              break;
            case Ar.BinaryExpression:
              this.recurse(e.left, void 0, void 0, (function(e) {
                  a = e
                })), this.recurse(e.right, void 0, void 0, (function(e) {
                  s = e
                })), l = "+" === e.operator ? this.plus(a, s) : "-" === e.operator ? this.ifDefined(a, 0) + e.operator + this.ifDefined(s, 0) : "(" +
                a + ")" + e.operator + "(" + s + ")", this.assign(t, l), r(l);
              break;
            case Ar.LogicalExpression:
              t = t || this.nextId(), d.recurse(e.left, t), d.if_("&&" === e.operator ? t : d.not(t), d.lazyRecurse(e.right, t)), r(t);
              break;
            case Ar.ConditionalExpression:
              t = t || this.nextId(), d.recurse(e.test, t), d.if_(t, d.lazyRecurse(e.alternate, t), d.lazyRecurse(e.consequent, t)), r(t);
              break;
            case Ar.Identifier:
              t = t || this.nextId(), n && (n.context = "inputs" === d.stage ? "s" : this.assign(this.nextId(), this.getHasOwnProperty("l", e.name) +
                "?l:s"), n.computed = !1, n.name = e.name), d.if_("inputs" === d.stage || d.not(d.getHasOwnProperty("l", e.name)), (function() {
                d.if_("inputs" === d.stage || "s", (function() {
                  i && 1 !== i && d.if_(d.isNull(d.nonComputedMember("s", e.name)), d.lazyAssign(d.nonComputedMember("s", e.name), "{}")),
                    d.assign(t, d.nonComputedMember("s", e.name))
                }))
              }), t && d.lazyAssign(t, d.nonComputedMember("l", e.name))), r(t);
              break;
            case Ar.MemberExpression:
              a = n && (n.context = this.nextId()) || this.nextId(), t = t || this.nextId(), d.recurse(e.object, a, void 0, (function() {
                d.if_(d.notNull(a), (function() {
                  e.computed ? (s = d.nextId(), d.recurse(e.property, s), d.getStringValue(s), i && 1 !== i && d.if_(d.not(d
                      .computedMember(a, s)), d.lazyAssign(d.computedMember(a, s), "{}")), l = d.computedMember(a, s), d.assign(t, l),
                    n && (n.computed = !0, n.name = s)) : (i && 1 !== i && d.if_(d.isNull(d.nonComputedMember(a, e.property.name)), d
                      .lazyAssign(d.nonComputedMember(a, e.property.name), "{}")), l = d.nonComputedMember(a, e.property.name), d
                    .assign(t, l), n && (n.computed = !1, n.name = e.property.name))
                }), (function() {
                  d.assign(t, "undefined")
                })), r(t)
              }), !!i);
              break;
            case Ar.CallExpression:
              t = t || this.nextId(), e.filter ? (s = d.filter(e.callee.name), u = [], C(e.arguments, (function(e) {
                var t = d.nextId();
                d.recurse(e, t), u.push(t)
              })), l = s + "(" + u.join(",") + ")", d.assign(t, l), r(t)) : (s = d.nextId(), a = {}, u = [], d.recurse(e.callee, s, a, (function() {
                d.if_(d.notNull(s), (function() {
                  C(e.arguments, (function(t) {
                      d.recurse(t, e.constant ? void 0 : d.nextId(), void 0, (function(e) {
                        u.push(e)
                      }))
                    })), l = a.name ? d.member(a.context, a.name, a.computed) + "(" + u.join(",") + ")" : s + "(" + u.join(",") + ")", d
                    .assign(t, l)
                }), (function() {
                  d.assign(t, "undefined")
                })), r(t)
              })));
              break;
            case Ar.AssignmentExpression:
              s = this.nextId(), a = {}, this.recurse(e.left, void 0, a, (function() {
                d.if_(d.notNull(a.context), (function() {
                  d.recurse(e.right, s), l = d.member(a.context, a.name, a.computed) + e.operator + s, d.assign(t, l), r(t || l)
                }))
              }), 1);
              break;
            case Ar.ArrayExpression:
              u = [], C(e.elements, (function(t) {
                d.recurse(t, e.constant ? void 0 : d.nextId(), void 0, (function(e) {
                  u.push(e)
                }))
              })), l = "[" + u.join(",") + "]", this.assign(t, l), r(t || l);
              break;
            case Ar.ObjectExpression:
              u = [], c = !1, C(e.properties, (function(e) {
                e.computed && (c = !0)
              })), c ? (t = t || this.nextId(), this.assign(t, "{}"), C(e.properties, (function(e) {
                e.computed ? (a = d.nextId(), d.recurse(e.key, a)) : a = e.key.type === Ar.Identifier ? e.key.name : "" + e.key.value, s = d
                  .nextId(), d.recurse(e.value, s), d.assign(d.member(t, a, e.computed), s)
              }))) : (C(e.properties, (function(t) {
                d.recurse(t.value, e.constant ? void 0 : d.nextId(), void 0, (function(e) {
                  u.push(d.escape(t.key.type === Ar.Identifier ? t.key.name : "" + t.key.value) + ":" + e)
                }))
              })), l = "{" + u.join(",") + "}", this.assign(t, l)), r(t || l);
              break;
            case Ar.ThisExpression:
              this.assign(t, "s"), r(t || "s");
              break;
            case Ar.LocalsExpression:
              this.assign(t, "l"), r(t || "l");
              break;
            case Ar.NGValueParameter:
              this.assign(t, "v"), r(t || "v")
          }
        },
        getHasOwnProperty: function(e, t) {
          var n = e + "." + t,
            r = this.current().own;
          return r.hasOwnProperty(n) || (r[n] = this.nextId(!1, e + "&&(" + this.escape(t) + " in " + e + ")")), r[n]
        },
        assign: function(e, t) {
          if (e) return this.current().body.push(e, "=", t, ";"), e
        },
        filter: function(e) {
          return this.state.filters.hasOwnProperty(e) || (this.state.filters[e] = this.nextId(!0)), this.state.filters[e]
        },
        ifDefined: function(e, t) {
          return "ifDefined(" + e + "," + this.escape(t) + ")"
        },
        plus: function(e, t) {
          return "plus(" + e + "," + t + ")"
        },
        return_: function(e) {
          this.current().body.push("return ", e, ";")
        },
        if_: function(e, t, n) {
          if (!0 === e) t();
          else {
            var r = this.current().body;
            r.push("if(", e, "){"), t(), r.push("}"), n && (r.push("else{"), n(), r.push("}"))
          }
        },
        not: function(e) {
          return "!(" + e + ")"
        },
        isNull: function(e) {
          return e + "==null"
        },
        notNull: function(e) {
          return e + "!=null"
        },
        nonComputedMember: function(e, t) {
          return /^[$_a-zA-Z][$_a-zA-Z0-9]*$/.test(t) ? e + "." + t : e + '["' + t.replace(/[^$_a-zA-Z0-9]/g, this.stringEscapeFn) + '"]'
        },
        computedMember: function(e, t) {
          return e + "[" + t + "]"
        },
        member: function(e, t, n) {
          return n ? this.computedMember(e, t) : this.nonComputedMember(e, t)
        },
        getStringValue: function(e) {
          this.assign(e, "getStringValue(" + e + ")")
        },
        lazyRecurse: function(e, t, n, r, i, o) {
          var a = this;
          return function() {
            a.recurse(e, t, n, r, i, o)
          }
        },
        lazyAssign: function(e, t) {
          var n = this;
          return function() {
            n.assign(e, t)
          }
        },
        stringEscapeRegex: /[^ a-zA-Z0-9]/g,
        stringEscapeFn: function(e) {
          return "\\u" + ("0000" + e.charCodeAt(0).toString(16)).slice(-4)
        },
        escape: function(e) {
          if (H(e)) return "'" + e.replace(this.stringEscapeRegex, this.stringEscapeFn) + "'";
          if (F(e)) return e.toString();
          if (!0 === e) return "true";
          if (!1 === e) return "false";
          if (null === e) return "null";
          if (void 0 === e) return "undefined";
          throw br("esc", "IMPOSSIBLE")
        },
        nextId: function(e, t) {
          var n = "v" + this.state.nextId++;
          return e || this.current().vars.push(n + (t ? "=" + t : "")), n
        },
        current: function() {
          return this.state[this.state.computing]
        }
      }, Pr.prototype = {
        compile: function(e) {
          var t, n, r = this;
          Dr(e, r.$filter), (t = Or(e)) && (n = this.recurse(t));
          var i, o = Mr(e.body);
          o && (i = [], C(o, (function(e, t) {
            var n = r.recurse(e);
            n.isPure = e.isPure, e.input = n, i.push(n), e.watchId = t
          })));
          var a = [];
          C(e.body, (function(e) {
            a.push(r.recurse(e.expression))
          }));
          var s = 0 === e.body.length ? q : 1 === e.body.length ? a[0] : function(e, t) {
            var n;
            return C(a, (function(r) {
              n = r(e, t)
            })), n
          };
          return n && (s.assign = function(e, t, r) {
            return n(e, r, t)
          }), i && (s.inputs = i), s
        },
        recurse: function(e, t, n) {
          var r, i, o, a = this;
          if (e.input) return this.inputs(e.input, e.watchId);
          switch (e.type) {
            case Ar.Literal:
              return this.value(e.value, t);
            case Ar.UnaryExpression:
              return i = this.recurse(e.argument), this["unary" + e.operator](i, t);
            case Ar.BinaryExpression:
            case Ar.LogicalExpression:
              return r = this.recurse(e.left), i = this.recurse(e.right), this["binary" + e.operator](r, i, t);
            case Ar.ConditionalExpression:
              return this["ternary?:"](this.recurse(e.test), this.recurse(e.alternate), this.recurse(e.consequent), t);
            case Ar.Identifier:
              return a.identifier(e.name, t, n);
            case Ar.MemberExpression:
              return r = this.recurse(e.object, !1, !!n), e.computed || (i = e.property.name), e.computed && (i = this.recurse(e.property)), e
                .computed ? this.computedMember(r, i, t, n) : this.nonComputedMember(r, i, t, n);
            case Ar.CallExpression:
              return o = [], C(e.arguments, (function(e) {
                o.push(a.recurse(e))
              })), e.filter && (i = this.$filter(e.callee.name)), e.filter || (i = this.recurse(e.callee, !0)), e.filter ? function(e, n, r, a) {
                for (var s = [], u = 0; u < o.length; ++u) s.push(o[u](e, n, r, a));
                var l = i.apply(void 0, s, a);
                return t ? {
                  context: void 0,
                  name: void 0,
                  value: l
                } : l
              } : function(e, n, r, a) {
                var s, u = i(e, n, r, a);
                if (null != u.value) {
                  for (var l = [], c = 0; c < o.length; ++c) l.push(o[c](e, n, r, a));
                  s = u.value.apply(u.context, l)
                }
                return t ? {
                  value: s
                } : s
              };
            case Ar.AssignmentExpression:
              return r = this.recurse(e.left, !0, 1), i = this.recurse(e.right),
                function(e, n, o, a) {
                  var s = r(e, n, o, a),
                    u = i(e, n, o, a);
                  return s.context[s.name] = u, t ? {
                    value: u
                  } : u
                };
            case Ar.ArrayExpression:
              return o = [], C(e.elements, (function(e) {
                  o.push(a.recurse(e))
                })),
                function(e, n, r, i) {
                  for (var a = [], s = 0; s < o.length; ++s) a.push(o[s](e, n, r, i));
                  return t ? {
                    value: a
                  } : a
                };
            case Ar.ObjectExpression:
              return o = [], C(e.properties, (function(e) {
                  e.computed ? o.push({
                    key: a.recurse(e.key),
                    computed: !0,
                    value: a.recurse(e.value)
                  }) : o.push({
                    key: e.key.type === Ar.Identifier ? e.key.name : "" + e.key.value,
                    computed: !1,
                    value: a.recurse(e.value)
                  })
                })),
                function(e, n, r, i) {
                  for (var a = {}, s = 0; s < o.length; ++s) o[s].computed ? a[o[s].key(e, n, r, i)] = o[s].value(e, n, r, i) : a[o[s].key] = o[s]
                    .value(e, n, r, i);
                  return t ? {
                    value: a
                  } : a
                };
            case Ar.ThisExpression:
              return function(e) {
                return t ? {
                  value: e
                } : e
              };
            case Ar.LocalsExpression:
              return function(e, n) {
                return t ? {
                  value: n
                } : n
              };
            case Ar.NGValueParameter:
              return function(e, n, r) {
                return t ? {
                  value: r
                } : r
              }
          }
        },
        "unary+": function(e, t) {
          return function(n, r, i, o) {
            var a = e(n, r, i, o);
            return a = _(a) ? +a : 0, t ? {
              value: a
            } : a
          }
        },
        "unary-": function(e, t) {
          return function(n, r, i, o) {
            var a = e(n, r, i, o);
            return a = _(a) ? -a : -0, t ? {
              value: a
            } : a
          }
        },
        "unary!": function(e, t) {
          return function(n, r, i, o) {
            var a = !e(n, r, i, o);
            return t ? {
              value: a
            } : a
          }
        },
        "binary+": function(e, t, n) {
          return function(r, i, o, a) {
            var s = kr(e(r, i, o, a), t(r, i, o, a));
            return n ? {
              value: s
            } : s
          }
        },
        "binary-": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a),
              u = t(r, i, o, a),
              l = (_(s) ? s : 0) - (_(u) ? u : 0);
            return n ? {
              value: l
            } : l
          }
        },
        "binary*": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) * t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary/": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) / t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary%": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) % t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary===": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) === t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary!==": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) !== t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary==": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) == t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary!=": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) != t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary<": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) < t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary>": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) > t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary<=": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) <= t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary>=": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) >= t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary&&": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) && t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "binary||": function(e, t, n) {
          return function(r, i, o, a) {
            var s = e(r, i, o, a) || t(r, i, o, a);
            return n ? {
              value: s
            } : s
          }
        },
        "ternary?:": function(e, t, n, r) {
          return function(i, o, a, s) {
            var u = e(i, o, a, s) ? t(i, o, a, s) : n(i, o, a, s);
            return r ? {
              value: u
            } : u
          }
        },
        value: function(e, t) {
          return function() {
            return t ? {
              context: void 0,
              name: void 0,
              value: e
            } : e
          }
        },
        identifier: function(e, t, n) {
          return function(r, i, o, a) {
            var s = i && e in i ? i : r;
            n && 1 !== n && s && null == s[e] && (s[e] = {});
            var u = s ? s[e] : void 0;
            return t ? {
              context: s,
              name: e,
              value: u
            } : u
          }
        },
        computedMember: function(e, t, n, r) {
          return function(i, o, a, s) {
            var u, l, c = e(i, o, a, s);
            return null != c && (u = wr(u = t(i, o, a, s)), r && 1 !== r && c && !c[u] && (c[u] = {}), l = c[u]), n ? {
              context: c,
              name: u,
              value: l
            } : l
          }
        },
        nonComputedMember: function(e, t, n, r) {
          return function(i, o, a, s) {
            var u = e(i, o, a, s);
            r && 1 !== r && u && null == u[t] && (u[t] = {});
            var l = null != u ? u[t] : void 0;
            return n ? {
              context: u,
              name: t,
              value: l
            } : l
          }
        },
        inputs: function(e, t) {
          return function(n, r, i, o) {
            return o ? o[t] : e(n, r, i)
          }
        }
      }, qr.prototype = {
        constructor: qr,
        parse: function(e) {
          var t = this.getAst(e),
            n = this.astCompiler.compile(t.ast);
          return n.literal = function(e) {
            return 0 === e.body.length || 1 === e.body.length && (e.body[0].expression.type === Ar.Literal || e.body[0].expression.type === Ar
              .ArrayExpression || e.body[0].expression.type === Ar.ObjectExpression)
          }(t.ast), n.constant = function(e) {
            return e.constant
          }(t.ast), n.oneTime = t.oneTime, n
        },
        getAst: function(e) {
          var t = !1;
          return ":" === (e = e.trim()).charAt(0) && ":" === e.charAt(1) && (t = !0, e = e.substring(2)), {
            ast: this.ast.ast(e),
            oneTime: t
          }
        }
      };
      var Wr = i("$sce"),
        zr = {
          HTML: "html",
          CSS: "css",
          MEDIA_URL: "mediaUrl",
          URL: "url",
          RESOURCE_URL: "resourceUrl",
          JS: "js"
        },
        Yr = /_([a-z])/g;

      function Kr(e) {
        return e.replace(Yr, Je)
      }

      function Qr(e) {
        var t = [];
        return _(e) && C(e, (function(e) {
          t.push(function(e) {
            if ("self" === e) return e;
            if (H(e)) {
              if (e.indexOf("***") > -1) throw Wr("iwcard", "Illegal sequence *** in string matcher.  String: {0}", e);
              return e = te(e).replace(/\\\*\\\*/g, ".*").replace(/\\\*/g, "[^:/.?&;]*"), new RegExp("^" + e + "$")
            }
            if (Y(e)) return new RegExp("^" + e.source + "$");
            throw Wr("imatcher", 'Matchers may only be "self", string patterns or RegExp objects')
          }(e))
        })), t
      }

      function Xr() {
        this.SCE_CONTEXTS = zr;
        var t = ["self"],
          n = [];
        this.trustedResourceUrlList = function(e) {
          return arguments.length && (t = Qr(e)), t
        }, Object.defineProperty(this, "resourceUrlWhitelist", {
          get: function() {
            return this.trustedResourceUrlList
          },
          set: function(e) {
            this.trustedResourceUrlList = e
          }
        }), this.bannedResourceUrlList = function(e) {
          return arguments.length && (n = Qr(e)), n
        }, Object.defineProperty(this, "resourceUrlBlacklist", {
          get: function() {
            return this.bannedResourceUrlList
          },
          set: function(e) {
            this.bannedResourceUrlList = e
          }
        }), this.$get = ["$injector", "$$sanitizeUri", function(r, i) {
          var o = function(e) {
            throw Wr("unsafe", "Attempting to use an unsafe value in a safe context.")
          };

          function a(t, n) {
            return "self" === t ? pi(n, li) || function(t) {
              return pi(t, function() {
                if (e.document.baseURI) return e.document.baseURI;
                si || ((si = e.document.createElement("a")).href = ".", si = si.cloneNode(!1));
                return si.href
              }())
            }(n) : !!t.exec(n.href)
          }

          function s(e) {
            var t = function(e) {
              this.$$unwrapTrustedValue = function() {
                return e
              }
            };
            return e && (t.prototype = new e), t.prototype.valueOf = function() {
              return this.$$unwrapTrustedValue()
            }, t.prototype.toString = function() {
              return this.$$unwrapTrustedValue().toString()
            }, t
          }
          r.has("$sanitize") && (o = r.get("$sanitize"));
          var u = s(),
            l = {};
          return l[zr.HTML] = s(u), l[zr.CSS] = s(u), l[zr.MEDIA_URL] = s(u), l[zr.URL] = s(l[zr.MEDIA_URL]), l[zr.JS] = s(u), l[zr.RESOURCE_URL] = s(
            l[zr.URL]), {
            trustAs: function(e, t) {
              var n = l.hasOwnProperty(e) ? l[e] : null;
              if (!n) throw Wr("icontext", "Attempted to trust a value in invalid context. Context: {0}; Value: {1}", e, t);
              if (null === t || j(t) || "" === t) return t;
              if ("string" != typeof t) throw Wr("itype", "Attempted to trust a non-string value in a content requiring a string: Context: {0}",
                e);
              return new n(t)
            },
            getTrusted: function(e, r) {
              if (null === r || j(r) || "" === r) return r;
              var s = l.hasOwnProperty(e) ? l[e] : null;
              if (s && r instanceof s) return r.$$unwrapTrustedValue();
              if (z(r.$$unwrapTrustedValue) && (r = r.$$unwrapTrustedValue()), e === zr.MEDIA_URL || e === zr.URL) return i(r.toString(), e === zr
                .MEDIA_URL);
              if (e === zr.RESOURCE_URL) {
                if (function(e) {
                    var r, i, o = di(e.toString()),
                      s = !1;
                    for (r = 0, i = t.length; r < i; r++)
                      if (a(t[r], o)) {
                        s = !0;
                        break
                      } if (s)
                      for (r = 0, i = n.length; r < i; r++)
                        if (a(n[r], o)) {
                          s = !1;
                          break
                        } return s
                  }(r)) return r;
                throw Wr("insecurl", "Blocked loading resource from url not allowed by $sceDelegate policy.  URL: {0}", r.toString())
              }
              if (e === zr.HTML) return o(r);
              throw Wr("unsafe", "Attempting to use an unsafe value in a safe context.")
            },
            valueOf: function(e) {
              return e instanceof u ? e.$$unwrapTrustedValue() : e
            }
          }
        }]
      }

      function Jr() {
        var e = !0;
        this.enabled = function(t) {
          return arguments.length && (e = !!t), e
        }, this.$get = ["$parse", "$sceDelegate", function(t, n) {
          if (e && o < 8) throw Wr("iequirks",
            "Strict Contextual Escaping does not support Internet Explorer version < 11 in quirks mode.  You can fix this by adding the text <!doctype html> to the top of your HTML document.  See http://docs.angularjs.org/api/ng.$sce for more information."
            );
          var r = Fe(zr);
          r.isEnabled = function() {
            return e
          }, r.trustAs = n.trustAs, r.getTrusted = n.getTrusted, r.valueOf = n.valueOf, e || (r.trustAs = r.getTrusted = function(e, t) {
            return t
          }, r.valueOf = I), r.parseAs = function(e, n) {
            var i = t(n);
            return i.literal && i.constant ? i : t(n, (function(t) {
              return r.getTrusted(e, t)
            }))
          };
          var i = r.parseAs,
            a = r.getTrusted,
            s = r.trustAs;
          return C(zr, (function(e, t) {
            var n = p(t);
            r[Kr("parse_as_" + n)] = function(t) {
              return i(e, t)
            }, r[Kr("get_trusted_" + n)] = function(t) {
              return a(e, t)
            }, r[Kr("trust_as_" + n)] = function(t) {
              return s(e, t)
            }
          })), r
        }]
      }

      function Zr() {
        this.$get = ["$window", "$document", function(e, t) {
          var n = {},
            r = !(!(e.nw && e.nw.process) && e.chrome && (e.chrome.app && e.chrome.app.runtime || !e.chrome.app && e.chrome.runtime && e.chrome
              .runtime.id)) && e.history && e.history.pushState,
            i = O((/android (\d+)/.exec(p((e.navigator || {}).userAgent)) || [])[1]),
            a = /Boxee/i.test((e.navigator || {}).userAgent),
            s = t[0] || {},
            u = s.body && s.body.style,
            l = !1,
            c = !1;
          return u && (l = !(!("transition" in u) && !("webkitTransition" in u)), c = !(!("animation" in u) && !("webkitAnimation" in u))), {
            history: !(!r || i < 4 || a),
            hasEvent: function(e) {
              if ("input" === e && o) return !1;
              if (j(n[e])) {
                var t = s.createElement("div");
                n[e] = "on" + e in t
              }
              return n[e]
            },
            csp: le(),
            transitions: l,
            animations: c,
            android: i
          }
        }]
      }

      function ei() {
        this.$get = L((function(e) {
          return new ti(e)
        }))
      }

      function ti(e) {
        var t = this,
          n = {},
          r = [],
          i = t.ALL_TASKS_TYPE = "$$all$$",
          o = t.DEFAULT_TASK_TYPE = "$$default$$";

        function a() {
          var e = r.pop();
          return e && e.cb
        }

        function s(e) {
          for (var t = r.length - 1; t >= 0; --t) {
            var n = r[t];
            if (n.type === e) return r.splice(t, 1), n.cb
          }
        }
        t.completeTask = function(t, r) {
          r = r || o;
          try {
            t()
          } finally {
            ! function(e) {
              n[e = e || o] && (n[e]--, n[i]--)
            }(r);
            var u = n[r],
              l = n[i];
            if (!l || !u)
              for (var c, d = l ? s : a; c = d(r);) try {
                c()
              } catch (t) {
                e.error(t)
              }
          }
        }, t.incTaskCount = function(e) {
          n[e = e || o] = (n[e] || 0) + 1, n[i] = (n[i] || 0) + 1
        }, t.notifyWhenNoPendingTasks = function(e, t) {
          n[t = t || i] ? r.push({
            type: t,
            cb: e
          }) : e()
        }
      }
      var ni = i("$templateRequest");

      function ri() {
        var e;
        this.httpOptions = function(t) {
          return t ? (e = t, this) : e
        }, this.$get = ["$exceptionHandler", "$templateCache", "$http", "$q", "$sce", function(t, n, r, i, o) {
          function a(s, u) {
            a.totalPendingRequests++, H(s) && !j(n.get(s)) || (s = o.getTrustedResourceUrl(s));
            var l = r.defaults && r.defaults.transformResponse;
            return G(l) ? l = l.filter((function(e) {
              return e !== Un
            })) : l === Un && (l = null), r.get(s, M({
              cache: n,
              transformResponse: l
            }, e)).finally((function() {
              a.totalPendingRequests--
            })).then((function(e) {
              return n.put(s, e.data)
            }), (function(e) {
              u || (e = ni("tpload", "Failed to load template: {0} (HTTP status: {1} {2})", s, e.status, e.statusText), t(e));
              return i.reject(e)
            }))
          }
          return a.totalPendingRequests = 0, a
        }]
      }

      function ii() {
        this.$get = ["$rootScope", "$browser", "$location", function(e, t, n) {
          var r = {
            findBindings: function(e, t, n) {
              var r = e.getElementsByClassName("ng-binding"),
                i = [];
              return C(r, (function(e) {
                var r = y.element(e).data("$binding");
                r && C(r, (function(r) {
                  n ? new RegExp("(^|\\s)" + te(t) + "(\\s|\\||$)").test(r) && i.push(e) : -1 !== r.indexOf(t) && i.push(e)
                }))
              })), i
            },
            findModels: function(e, t, n) {
              for (var r = ["ng-", "data-ng-", "ng\\:"], i = 0; i < r.length; ++i) {
                var o = "[" + r[i] + "model" + (n ? "=" : "*=") + '"' + t + '"]',
                  a = e.querySelectorAll(o);
                if (a.length) return a
              }
            },
            getLocation: function() {
              return n.url()
            },
            setLocation: function(t) {
              t !== n.url() && (n.url(t), e.$digest())
            },
            whenStable: function(e) {
              t.notifyWhenNoOutstandingRequests(e)
            }
          };
          return r
        }]
      }
      var oi = i("$timeout");

      function ai() {
        this.$get = ["$rootScope", "$browser", "$q", "$$q", "$exceptionHandler", function(e, t, n, r, i) {
          var o = {};

          function a(a, s, u) {
            z(a) || (u = s, s = a, a = q);
            var l, c = pe(arguments, 3),
              d = _(u) && !u,
              p = (d ? r : n).defer(),
              f = p.promise;
            return l = t.defer((function() {
              try {
                p.resolve(a.apply(null, c))
              } catch (e) {
                p.reject(e), i(e)
              } finally {
                delete o[f.$$timeoutId]
              }
              d || e.$apply()
            }), s, "$timeout"), f.$$timeoutId = l, o[l] = p, f
          }
          return a.cancel = function(e) {
            if (!e) return !1;
            if (!e.hasOwnProperty("$$timeoutId")) throw oi("badprom",
              "`$timeout.cancel()` called with a promise that was not generated by `$timeout()`.");
            if (!o.hasOwnProperty(e.$$timeoutId)) return !1;
            var n = e.$$timeoutId,
              r = o[n];
            return Hr(r.promise), r.reject("canceled"), delete o[n], t.defer.cancel(n)
          }, a
        }]
      }
      var si, ui = e.document.createElement("a"),
        li = di(e.location.href);
      ui.href = "http://[::1]";
      var ci = "[::1]" === ui.hostname;

      function di(e) {
        if (!H(e)) return e;
        var t = e;
        o && (ui.setAttribute("href", t), t = ui.href), ui.setAttribute("href", t);
        var n = ui.hostname;
        return !ci && n.indexOf(":") > -1 && (n = "[" + n + "]"), {
          href: ui.href,
          protocol: ui.protocol ? ui.protocol.replace(/:$/, "") : "",
          host: ui.host,
          search: ui.search ? ui.search.replace(/^\?/, "") : "",
          hash: ui.hash ? ui.hash.replace(/^#/, "") : "",
          hostname: n,
          port: ui.port,
          pathname: "/" === ui.pathname.charAt(0) ? ui.pathname : "/" + ui.pathname
        }
      }

      function pi(e, t) {
        return e = di(e), t = di(t), e.protocol === t.protocol && e.host === t.host
      }

      function fi() {
        this.$get = L(e)
      }

      function hi(e) {
        var t = e[0] || {},
          n = {},
          r = "";

        function i(e) {
          try {
            return decodeURIComponent(e)
          } catch (t) {
            return e
          }
        }
        return function() {
          var e, o, a, s, u, l = function(e) {
            try {
              return e.cookie || ""
            } catch (e) {
              return ""
            }
          }(t);
          if (l !== r)
            for (e = (r = l).split("; "), n = {}, a = 0; a < e.length; a++)(s = (o = e[a]).indexOf("=")) > 0 && (u = i(o.substring(0, s)), j(n[u]) && (
              n[u] = i(o.substring(s + 1))));
          return n
        }
      }

      function gi() {
        this.$get = hi
      }

      function mi(e) {
        var t = "Filter";

        function n(r, i) {
          if (V(r)) {
            var o = {};
            return C(r, (function(e, t) {
              o[t] = n(t, e)
            })), o
          }
          return e.factory(r + t, i)
        }
        this.register = n, this.$get = ["$injector", function(e) {
            return function(n) {
              return e.get(n + t)
            }
          }], n("currency", yi), n("date", Ni), n("filter", vi), n("json", Pi), n("limitTo", Li), n("lowercase", qi), n("number", wi), n("orderBy", ji),
          n("uppercase", Ii)
      }

      function vi() {
        return function(e, t, n, r) {
          if (!x(e)) {
            if (null == e) return e;
            throw i("filter")("notarray", "Expected array but received: {0}", e)
          }
          var o, a;
          switch (r = r || "$", bi(t)) {
            case "function":
              o = t;
              break;
            case "boolean":
            case "null":
            case "number":
            case "string":
              a = !0;
            case "object":
              o = function(e, t, n, r) {
                var i = V(e) && n in e;
                !0 === t ? t = ue : z(t) || (t = function(e, t) {
                  return !j(e) && (null === e || null === t ? e === t : !(V(t) || V(e) && !R(e)) && (e = p("" + e), t = p("" + t), -1 !== e
                    .indexOf(t)))
                });
                return function(o) {
                  return i && !V(o) ? $i(o, e[n], t, n, !1) : $i(o, e, t, n, r)
                }
              }(t, n, r, a);
              break;
            default:
              return e
          }
          return Array.prototype.filter.call(e, o)
        }
      }

      function $i(e, t, n, r, i, o) {
        var a = bi(e),
          s = bi(t);
        if ("string" === s && "!" === t.charAt(0)) return !$i(e, t.substring(1), n, r, i);
        if (G(e)) return e.some((function(e) {
          return $i(e, t, n, r, i)
        }));
        switch (a) {
          case "object":
            var u;
            if (i) {
              for (u in e)
                if (u.charAt && "$" !== u.charAt(0) && $i(e[u], t, n, r, !0)) return !0;
              return !o && $i(e, t, n, r, !1)
            }
            if ("object" === s) {
              for (u in t) {
                var l = t[u];
                if (!z(l) && !j(l)) {
                  var c = u === r;
                  if (!$i(c ? e : e[u], l, n, r, c, c)) return !1
                }
              }
              return !0
            }
            return n(e, t);
          case "function":
            return !1;
          default:
            return n(e, t)
        }
      }

      function bi(e) {
        return null === e ? "null" : typeof e
      }
      hi.$inject = ["$document"], mi.$inject = ["$provide"];

      function yi(e) {
        var t = e.NUMBER_FORMATS;
        return function(e, n, r) {
          j(n) && (n = t.CURRENCY_SYM), j(r) && (r = t.PATTERNS[1].maxFrac);
          var i = n ? /\u00A4/g : /\s*\u00A4\s*/g;
          return null == e ? e : xi(e, t.PATTERNS[1], t.GROUP_SEP, t.DECIMAL_SEP, r).replace(i, n)
        }
      }

      function wi(e) {
        var t = e.NUMBER_FORMATS;
        return function(e, n) {
          return null == e ? e : xi(e, t.PATTERNS[0], t.GROUP_SEP, t.DECIMAL_SEP, n)
        }
      }

      function xi(e, t, n, r, i) {
        if (!H(e) && !F(e) || isNaN(e)) return "";
        var o, a = !isFinite(e),
          s = !1,
          u = Math.abs(e) + "",
          l = "";
        if (a) l = "∞";
        else {
          (function(e, t, n, r) {
            var i = e.d,
              o = i.length - e.i,
              a = (t = j(t) ? Math.min(Math.max(n, o), r) : +t) + e.i,
              s = i[a];
            if (a > 0) {
              i.splice(Math.max(e.i, a));
              for (var u = a; u < i.length; u++) i[u] = 0
            } else {
              o = Math.max(0, o), e.i = 1, i.length = Math.max(1, a = t + 1), i[0] = 0;
              for (var l = 1; l < a; l++) i[l] = 0
            }
            if (s >= 5)
              if (a - 1 < 0) {
                for (var c = 0; c > a; c--) i.unshift(0), e.i++;
                i.unshift(1), e.i++
              } else i[a - 1]++;
            for (; o < Math.max(0, t); o++) i.push(0);
            var d = i.reduceRight((function(e, t, n, r) {
              return t += e, r[n] = t % 10, Math.floor(t / 10)
            }), 0);
            d && (i.unshift(d), e.i++)
          })(o = function(e) {
            var t, n, r, i, o, a = 0;
            for ((n = e.indexOf(".")) > -1 && (e = e.replace(".", "")), (r = e.search(/e/i)) > 0 ? (n < 0 && (n = r), n += +e.slice(r + 1), e = e
                .substring(0, r)) : n < 0 && (n = e.length), r = 0;
              "0" === e.charAt(r); r++);
            if (r === (o = e.length)) t = [0], n = 1;
            else {
              for (o--;
                "0" === e.charAt(o);) o--;
              for (n -= r, t = [], i = 0; r <= o; r++, i++) t[i] = +e.charAt(r)
            }
            return n > 22 && (t = t.splice(0, 21), a = n - 1, n = 1), {
              d: t,
              e: a,
              i: n
            }
          }(u), i, t.minFrac, t.maxFrac);
          var c = o.d,
            d = o.i,
            p = o.e,
            f = [];
          for (s = c.reduce((function(e, t) {
              return e && !t
            }), !0); d < 0;) c.unshift(0), d++;
          d > 0 ? f = c.splice(d, c.length) : (f = c, c = [0]);
          var h = [];
          for (c.length >= t.lgSize && h.unshift(c.splice(-t.lgSize, c.length).join("")); c.length > t.gSize;) h.unshift(c.splice(-t.gSize, c.length)
            .join(""));
          c.length && h.unshift(c.join("")), l = h.join(n), f.length && (l += r + f.join("")), p && (l += "e+" + p)
        }
        return e < 0 && !s ? t.negPre + l + t.negSuf : t.posPre + l + t.posSuf
      }

      function Ci(e, t, n, r) {
        var i = "";
        for ((e < 0 || r && e <= 0) && (r ? e = 1 - e : (e = -e, i = "-")), e = "" + e; e.length < t;) e = "0" + e;
        return n && (e = e.substr(e.length - t)), i + e
      }

      function Si(e, t, n, r, i) {
        return n = n || 0,
          function(o) {
            var a = o["get" + e]();
            return (n > 0 || a > -n) && (a += n), 0 === a && -12 === n && (a = 12), Ci(a, t, r, i)
          }
      }

      function Ai(e, t, n) {
        return function(r, i) {
          var o = r["get" + e]();
          return i[f((n ? "STANDALONE" : "") + (t ? "SHORT" : "") + e)][o]
        }
      }

      function Ti(e) {
        var t = new Date(e, 0, 1).getDay();
        return new Date(e, 0, (t <= 4 ? 5 : 12) - t)
      }

      function ki(e) {
        return function(t) {
          var n, r = Ti(t.getFullYear()),
            i = +(n = t, new Date(n.getFullYear(), n.getMonth(), n.getDate() + (4 - n.getDay()))) - +r;
          return Ci(1 + Math.round(i / 6048e5), e)
        }
      }

      function Di(e, t) {
        return e.getFullYear() <= 0 ? t.ERAS[0] : t.ERAS[1]
      }
      yi.$inject = ["$locale"], wi.$inject = ["$locale"];
      var Mi = {
          yyyy: Si("FullYear", 4, 0, !1, !0),
          yy: Si("FullYear", 2, 0, !0, !0),
          y: Si("FullYear", 1, 0, !1, !0),
          MMMM: Ai("Month"),
          MMM: Ai("Month", !0),
          MM: Si("Month", 2, 1),
          M: Si("Month", 1, 1),
          LLLL: Ai("Month", !1, !0),
          dd: Si("Date", 2),
          d: Si("Date", 1),
          HH: Si("Hours", 2),
          H: Si("Hours", 1),
          hh: Si("Hours", 2, -12),
          h: Si("Hours", 1, -12),
          mm: Si("Minutes", 2),
          m: Si("Minutes", 1),
          ss: Si("Seconds", 2),
          s: Si("Seconds", 1),
          sss: Si("Milliseconds", 3),
          EEEE: Ai("Day"),
          EEE: Ai("Day", !0),
          a: function(e, t) {
            return e.getHours() < 12 ? t.AMPMS[0] : t.AMPMS[1]
          },
          Z: function(e, t, n) {
            var r = -1 * n,
              i = r >= 0 ? "+" : "";
            return i += Ci(Math[r > 0 ? "floor" : "ceil"](r / 60), 2) + Ci(Math.abs(r % 60), 2)
          },
          ww: ki(2),
          w: ki(1),
          G: Di,
          GG: Di,
          GGG: Di,
          GGGG: function(e, t) {
            return e.getFullYear() <= 0 ? t.ERANAMES[0] : t.ERANAMES[1]
          }
        },
        Ei = /((?:[^yMLdHhmsaZEwG']+)|(?:'(?:[^']|'')*')|(?:E+|y+|M+|L+|d+|H+|h+|m+|s+|a|Z|G+|w+))([\s\S]*)/,
        Oi = /^-?\d+$/;

      function Ni(e) {
        var t = /^(\d{4})-?(\d\d)-?(\d\d)(?:T(\d\d)(?::?(\d\d)(?::?(\d\d)(?:\.(\d+))?)?)?(Z|([+-])(\d\d):?(\d\d))?)?$/;
        return function(n, r, i) {
          var o, a, s = "",
            u = [];
          if (r = r || "mediumDate", r = e.DATETIME_FORMATS[r] || r, H(n) && (n = Oi.test(n) ? O(n) : function(e) {
              var n;
              if (n = e.match(t)) {
                var r = new Date(0),
                  i = 0,
                  o = 0,
                  a = n[8] ? r.setUTCFullYear : r.setFullYear,
                  s = n[8] ? r.setUTCHours : r.setHours;
                n[9] && (i = O(n[9] + n[10]), o = O(n[9] + n[11])), a.call(r, O(n[1]), O(n[2]) - 1, O(n[3]));
                var u = O(n[4] || 0) - i,
                  l = O(n[5] || 0) - o,
                  c = O(n[6] || 0),
                  d = Math.round(1e3 * parseFloat("0." + (n[7] || 0)));
                return s.call(r, u, l, c, d), r
              }
              return e
            }(n)), F(n) && (n = new Date(n)), !B(n) || !isFinite(n.getTime())) return n;
          for (; r;)(a = Ei.exec(r)) ? r = (u = de(u, a, 1)).pop() : (u.push(r), r = null);
          var l = n.getTimezoneOffset();
          return i && (l = $e(i, l), n = ye(n, i, !0)), C(u, (function(t) {
            s += (o = Mi[t]) ? o(n, e.DATETIME_FORMATS, l) : "''" === t ? "'" : t.replace(/(^'|'$)/g, "").replace(/''/g, "'")
          })), s
        }
      }

      function Pi() {
        return function(e, t) {
          return j(t) && (t = 2), ge(e, t)
        }
      }
      Ni.$inject = ["$locale"];
      var qi = L(p),
        Ii = L(f);

      function Li() {
        return function(e, t, n) {
          return t = Math.abs(Number(t)) === 1 / 0 ? Number(t) : O(t), N(t) ? e : (F(e) && (e = e.toString()), x(e) ? (n = (n = !n || isNaN(n) ? 0 : O(
            n)) < 0 ? Math.max(0, e.length + n) : n, t >= 0 ? Ri(e, n, n + t) : 0 === n ? Ri(e, t, e.length) : Ri(e, Math.max(0, n + t), n)) : e)
        }
      }

      function Ri(e, t, n) {
        return H(e) ? e.slice(t, n) : h.call(e, t, n)
      }

      function ji(e) {
        return function(r, o, a, s) {
          if (null == r) return r;
          if (!x(r)) throw i("orderBy")("notarray", "Expected array but received: {0}", r);
          G(o) || (o = [o]), 0 === o.length && (o = ["+"]);
          var u = o.map((function(t) {
              var n = 1,
                r = I;
              if (z(t)) r = t;
              else if (H(t) && ("+" !== t.charAt(0) && "-" !== t.charAt(0) || (n = "-" === t.charAt(0) ? -1 : 1, t = t.substring(1)), "" !== t && (
                  r = e(t)).constant)) {
                var i = r();
                r = function(e) {
                  return e[i]
                }
              }
              return {
                get: r,
                descending: n
              }
            })),
            l = a ? -1 : 1,
            c = z(s) ? s : n,
            d = Array.prototype.map.call(r, (function(e, n) {
              return {
                value: e,
                tieBreaker: {
                  value: n,
                  type: "number",
                  index: n
                },
                predicateValues: u.map((function(r) {
                  return function(e, n) {
                    var r = typeof e;
                    null === e ? r = "null" : "object" === r && (e = function(e) {
                      if (z(e.valueOf) && t(e = e.valueOf())) return e;
                      if (R(e) && t(e = e.toString())) return e;
                      return e
                    }(e));
                    return {
                      value: e,
                      type: r,
                      index: n
                    }
                  }(r.get(e), n)
                }))
              }
            }));
          return d.sort((function(e, t) {
            for (var r = 0, i = u.length; r < i; r++) {
              var o = c(e.predicateValues[r], t.predicateValues[r]);
              if (o) return o * u[r].descending * l
            }
            return (c(e.tieBreaker, t.tieBreaker) || n(e.tieBreaker, t.tieBreaker)) * l
          })), r = d.map((function(e) {
            return e.value
          }))
        };

        function t(e) {
          switch (typeof e) {
            case "number":
            case "boolean":
            case "string":
              return !0;
            default:
              return !1
          }
        }

        function n(e, t) {
          var n = 0,
            r = e.type,
            i = t.type;
          if (r === i) {
            var o = e.value,
              a = t.value;
            "string" === r ? (o = o.toLowerCase(), a = a.toLowerCase()) : "object" === r && (V(o) && (o = e.index), V(a) && (a = t.index)), o !== a && (
              n = o < a ? -1 : 1)
          } else n = "undefined" === r ? 1 : "undefined" === i ? -1 : "null" === r ? 1 : "null" === i || r < i ? -1 : 1;
          return n
        }
      }

      function _i(e) {
        return z(e) && (e = {
          link: e
        }), e.restrict = e.restrict || "AC", L(e)
      }
      ji.$inject = ["$parse"];
      var Vi = L({
          restrict: "E",
          compile: function(e, t) {
            if (!t.href && !t.xlinkHref) return function(e, t) {
              if ("a" === t[0].nodeName.toLowerCase()) {
                var n = "[object SVGAnimatedString]" === v.call(t.prop("href")) ? "xlink:href" : "href";
                t.on("click", (function(e) {
                  t.attr(n) || e.preventDefault()
                }))
              }
            }
          }
        }),
        Ui = {};
      C(Nt, (function(e, t) {
        if ("multiple" !== e) {
          var n = wn("ng-" + t),
            r = i;
          "checked" === e && (r = function(e, t, r) {
            r.ngModel !== r[n] && i(e, 0, r)
          }), Ui[n] = function() {
            return {
              restrict: "A",
              priority: 100,
              link: r
            }
          }
        }

        function i(e, r, i) {
          e.$watch(i[n], (function(e) {
            i.$set(t, !!e)
          }))
        }
      })), C(qt, (function(e, t) {
        Ui[t] = function() {
          return {
            priority: 100,
            link: function(e, n, r) {
              if ("ngPattern" === t && "/" === r.ngPattern.charAt(0)) {
                var i = r.ngPattern.match(l);
                if (i) return void r.$set("ngPattern", new RegExp(i[1], i[2]))
              }
              e.$watch(r[t], (function(e) {
                r.$set(t, e)
              }))
            }
          }
        }
      })), C(["src", "srcset", "href"], (function(e) {
        var t = wn("ng-" + e);
        Ui[t] = ["$sce", function(n) {
          return {
            priority: 99,
            link: function(r, i, a) {
              var s = e,
                u = e;
              "href" === e && "[object SVGAnimatedString]" === v.call(i.prop("href")) && (u = "xlinkHref", a.$attr[u] = "xlink:href", s =
                null), a.$set(t, n.getTrustedMediaUrl(a[t])), a.$observe(t, (function(t) {
                t ? (a.$set(u, t), o && s && i.prop(s, a[u])) : "href" === e && a.$set(u, null)
              }))
            }
          }
        }]
      }));
      var Hi = {
          $addControl: q,
          $getControls: L([]),
          $$renameControl: function(e, t) {
            e.$name = t
          },
          $removeControl: q,
          $setValidity: q,
          $setDirty: q,
          $setPristine: q,
          $setSubmitted: q,
          $$setSubmitted: q
        },
        Fi = "ng-pending",
        Bi = "ng-submitted";

      function Gi(e, t, n, r, i) {
        this.$$controls = [], this.$error = {}, this.$$success = {}, this.$pending = void 0, this.$name = i(t.name || t.ngForm || "")(n), this.$dirty = !
          1, this.$pristine = !0, this.$valid = !0, this.$invalid = !1, this.$submitted = !1, this.$$parentForm = Hi, this.$$element = e, this.$$animate =
          r, Ki(this)
      }
      Gi.$inject = ["$element", "$attrs", "$scope", "$animate", "$interpolate"], Gi.prototype = {
        $rollbackViewValue: function() {
          C(this.$$controls, (function(e) {
            e.$rollbackViewValue()
          }))
        },
        $commitViewValue: function() {
          C(this.$$controls, (function(e) {
            e.$commitViewValue()
          }))
        },
        $addControl: function(e) {
          je(e.$name, "input"), this.$$controls.push(e), e.$name && (this[e.$name] = e), e.$$parentForm = this
        },
        $getControls: function() {
          return Fe(this.$$controls)
        },
        $$renameControl: function(e, t) {
          var n = e.$name;
          this[n] === e && delete this[n], this[t] = e, e.$name = t
        },
        $removeControl: function(e) {
          e.$name && this[e.$name] === e && delete this[e.$name], C(this.$pending, (function(t, n) {
            this.$setValidity(n, null, e)
          }), this), C(this.$error, (function(t, n) {
            this.$setValidity(n, null, e)
          }), this), C(this.$$success, (function(t, n) {
            this.$setValidity(n, null, e)
          }), this), oe(this.$$controls, e), e.$$parentForm = Hi
        },
        $setDirty: function() {
          this.$$animate.removeClass(this.$$element, Go), this.$$animate.addClass(this.$$element, Wo), this.$dirty = !0, this.$pristine = !1, this
            .$$parentForm.$setDirty()
        },
        $setPristine: function() {
          this.$$animate.setClass(this.$$element, Go, Wo + " " + Bi), this.$dirty = !1, this.$pristine = !0, this.$submitted = !1, C(this.$$controls,
            (function(e) {
              e.$setPristine()
            }))
        },
        $setUntouched: function() {
          C(this.$$controls, (function(e) {
            e.$setUntouched()
          }))
        },
        $setSubmitted: function() {
          for (var e = this; e.$$parentForm && e.$$parentForm !== Hi;) e = e.$$parentForm;
          e.$$setSubmitted()
        },
        $$setSubmitted: function() {
          this.$$animate.addClass(this.$$element, Bi), this.$submitted = !0, C(this.$$controls, (function(e) {
            e.$$setSubmitted && e.$$setSubmitted()
          }))
        }
      }, Qi({
        clazz: Gi,
        set: function(e, t, n) {
          var r = e[t];
          r ? -1 === r.indexOf(n) && r.push(n) : e[t] = [n]
        },
        unset: function(e, t, n) {
          var r = e[t];
          r && (oe(r, n), 0 === r.length && delete e[t])
        }
      });
      var Wi = function(e) {
          return ["$timeout", "$parse", function(t, n) {
            return {
              name: "form",
              restrict: e ? "EAC" : "E",
              require: ["form", "^^?form"],
              controller: Gi,
              compile: function(n, i) {
                n.addClass(Go).addClass(Fo);
                var o = i.name ? "name" : !(!e || !i.ngForm) && "ngForm";
                return {
                  pre: function(e, n, i, a) {
                    var s = a[0];
                    if (!("action" in i)) {
                      var u = function(t) {
                        e.$apply((function() {
                          s.$commitViewValue(), s.$setSubmitted()
                        })), t.preventDefault()
                      };
                      n[0].addEventListener("submit", u), n.on("$destroy", (function() {
                        t((function() {
                          n[0].removeEventListener("submit", u)
                        }), 0, !1)
                      }))
                    }(a[1] || s.$$parentForm).$addControl(s);
                    var l = o ? r(s.$name) : q;
                    o && (l(e, s), i.$observe(o, (function(t) {
                      s.$name !== t && (l(e, void 0), s.$$parentForm.$$renameControl(s, t), (l = r(s.$name))(e, s))
                    }))), n.on("$destroy", (function() {
                      s.$$parentForm.$removeControl(s), l(e, void 0), M(s, Hi)
                    }))
                  }
                }
              }
            };

            function r(e) {
              return "" === e ? n('this[""]').assign : n(e).assign || q
            }
          }]
        },
        zi = Wi(),
        Yi = Wi(!0);

      function Ki(e) {
        e.$$classCache = {}, e.$$classCache[Bo] = !(e.$$classCache[Fo] = e.$$element.hasClass(Fo))
      }

      function Qi(e) {
        var t = e.clazz,
          n = e.set,
          r = e.unset;

        function i(e, t, n) {
          n && !e.$$classCache[t] ? (e.$$animate.addClass(e.$$element, t), e.$$classCache[t] = !0) : !n && e.$$classCache[t] && (e.$$animate.removeClass(e
            .$$element, t), e.$$classCache[t] = !1)
        }

        function o(e, t, n) {
          t = t ? "-" + Pe(t, "-") : "", i(e, Fo + t, !0 === n), i(e, Bo + t, !1 === n)
        }
        t.prototype.$setValidity = function(e, t, a) {
          var s;
          j(t) ? function(e, t, r, i) {
              e[t] || (e[t] = {});
              n(e[t], r, i)
            }(this, "$pending", e, a) : function(e, t, n, i) {
              e[t] && r(e[t], n, i);
              Xi(e[t]) && (e[t] = void 0)
            }(this, "$pending", e, a), X(t) ? t ? (r(this.$error, e, a), n(this.$$success, e, a)) : (n(this.$error, e, a), r(this.$$success, e, a)) : (
              r(this.$error, e, a), r(this.$$success, e, a)), this.$pending ? (i(this, Fi, !0), this.$valid = this.$invalid = void 0, o(this, "",
            null)) : (i(this, Fi, !1), this.$valid = Xi(this.$error), this.$invalid = !this.$valid, o(this, "", this.$valid)), o(this, e, s = this
              .$pending && this.$pending[e] ? void 0 : !this.$error[e] && (!!this.$$success[e] || null)), this.$$parentForm.$setValidity(e, s, this)
        }
      }

      function Xi(e) {
        if (e)
          for (var t in e)
            if (e.hasOwnProperty(t)) return !1;
        return !0
      }
      var Ji = /^\d{4,}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d:[0-5]\d\.\d+(?:[+-][0-2]\d:[0-5]\d|Z)$/,
        Zi = /^[a-z][a-z\d.+-]*:\/*(?:[^:@]+(?::[^@]+)?@)?(?:[^\s:/?#]+|\[[a-f\d:]+])(?::\d+)?(?:\/[^?#]*)?(?:\?[^#]*)?(?:#.*)?$/i,
        eo =
        /^(?=.{1,254}$)(?=.{1,64}@)[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+(\.[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+)*@[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?(\.[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?)*$/,
        to = /^\s*(-|\+)?(\d+|(\d*(\.\d*)))([eE][+-]?\d+)?\s*$/,
        no = /^(\d{4,})-(\d{2})-(\d{2})$/,
        ro = /^(\d{4,})-(\d\d)-(\d\d)T(\d\d):(\d\d)(?::(\d\d)(\.\d{1,3})?)?$/,
        io = /^(\d{4,})-W(\d\d)$/,
        oo = /^(\d{4,})-(\d\d)$/,
        ao = /^(\d\d):(\d\d)(?::(\d\d)(\.\d{1,3})?)?$/,
        so = Ve();
      C("date,datetime-local,month,time,week".split(","), (function(e) {
        so[e] = !0
      }));
      var uo = {
        text: function(e, t, n, r, i, o) {
          co(e, t, n, r, i, o), lo(r)
        },
        date: fo("date", no, po(no, ["yyyy", "MM", "dd"]), "yyyy-MM-dd"),
        "datetime-local": fo("datetimelocal", ro, po(ro, ["yyyy", "MM", "dd", "HH", "mm", "ss", "sss"]), "yyyy-MM-ddTHH:mm:ss.sss"),
        time: fo("time", ao, po(ao, ["HH", "mm", "ss", "sss"]), "HH:mm:ss.sss"),
        week: fo("week", io, (function(e, t) {
          if (B(e)) return e;
          if (H(e)) {
            io.lastIndex = 0;
            var n = io.exec(e);
            if (n) {
              var r = +n[1],
                i = +n[2],
                o = 0,
                a = 0,
                s = 0,
                u = 0,
                l = Ti(r),
                c = 7 * (i - 1);
              return t && (o = t.getHours(), a = t.getMinutes(), s = t.getSeconds(), u = t.getMilliseconds()), new Date(r, 0, l.getDate() + c, o,
                a, s, u)
            }
          }
          return NaN
        }), "yyyy-Www"),
        month: fo("month", oo, po(oo, ["yyyy", "MM"]), "yyyy-MM"),
        number: function(e, t, n, r, i, o, a, s) {
          var u;
          if (ho(e, t, n, r, "number"), go(r), co(e, t, n, r, i, o), _(n.min) || n.ngMin) {
            var l = n.min || s(n.ngMin)(e);
            u = mo(l), r.$validators.min = function(e, t) {
              return r.$isEmpty(t) || j(u) || t >= u
            }, n.$observe("min", (function(e) {
              e !== l && (u = mo(e), l = e, r.$validate())
            }))
          }
          if (_(n.max) || n.ngMax) {
            var c = n.max || s(n.ngMax)(e),
              d = mo(c);
            r.$validators.max = function(e, t) {
              return r.$isEmpty(t) || j(d) || t <= d
            }, n.$observe("max", (function(e) {
              e !== c && (d = mo(e), c = e, r.$validate())
            }))
          }
          if (_(n.step) || n.ngStep) {
            var p = n.step || s(n.ngStep)(e),
              f = mo(p);
            r.$validators.step = function(e, t) {
              return r.$isEmpty(t) || j(f) || bo(t, u || 0, f)
            }, n.$observe("step", (function(e) {
              e !== p && (f = mo(e), p = e, r.$validate())
            }))
          }
        },
        url: function(e, t, n, r, i, o) {
          co(e, t, n, r, i, o), lo(r), r.$validators.url = function(e, t) {
            var n = e || t;
            return r.$isEmpty(n) || Zi.test(n)
          }
        },
        email: function(e, t, n, r, i, o) {
          co(e, t, n, r, i, o), lo(r), r.$validators.email = function(e, t) {
            var n = e || t;
            return r.$isEmpty(n) || eo.test(n)
          }
        },
        radio: function(e, t, n, r) {
          var i = !n.ngTrim || "false" !== ee(n.ngTrim);
          j(n.name) && t.attr("name", T());
          t.on("change", (function(e) {
            var o;
            t[0].checked && (o = n.value, i && (o = ee(o)), r.$setViewValue(o, e && e.type))
          })), r.$render = function() {
            var e = n.value;
            i && (e = ee(e)), t[0].checked = e === r.$viewValue
          }, n.$observe("value", r.$render)
        },
        range: function(e, t, n, r, i, o) {
          ho(e, t, n, r, "range"), go(r), co(e, t, n, r, i, o);
          var a = r.$$hasNativeValidators && "range" === t[0].type,
            s = a ? 0 : void 0,
            u = a ? 100 : void 0,
            l = a ? 1 : void 0,
            c = t[0].validity,
            d = _(n.min),
            p = _(n.max),
            f = _(n.step),
            h = r.$render;
          r.$render = a && _(c.rangeUnderflow) && _(c.rangeOverflow) ? function() {
            h(), r.$setViewValue(t.val())
          } : h, d && (s = mo(n.min), r.$validators.min = a ? function() {
            return !0
          } : function(e, t) {
            return r.$isEmpty(t) || j(s) || t >= s
          }, g("min", (function(e) {
            if (s = mo(e), N(r.$modelValue)) return;
            if (a) {
              var n = t.val();
              s > n && (n = s, t.val(n)), r.$setViewValue(n)
            } else r.$validate()
          })));
          p && (u = mo(n.max), r.$validators.max = a ? function() {
            return !0
          } : function(e, t) {
            return r.$isEmpty(t) || j(u) || t <= u
          }, g("max", (function(e) {
            if (u = mo(e), N(r.$modelValue)) return;
            if (a) {
              var n = t.val();
              u < n && (t.val(u), n = u < s ? s : u), r.$setViewValue(n)
            } else r.$validate()
          })));
          f && (l = mo(n.step), r.$validators.step = a ? function() {
            return !c.stepMismatch
          } : function(e, t) {
            return r.$isEmpty(t) || j(l) || bo(t, s || 0, l)
          }, g("step", (function(e) {
            if (l = mo(e), N(r.$modelValue)) return;
            a ? r.$viewValue !== t.val() && r.$setViewValue(t.val()) : r.$validate()
          })));

          function g(e, r) {
            t.attr(e, n[e]);
            var i = n[e];
            n.$observe(e, (function(e) {
              e !== i && (i = e, r(e))
            }))
          }
        },
        checkbox: function(e, t, n, r, i, o, a, s) {
          var u = yo(s, e, "ngTrueValue", n.ngTrueValue, !0),
            l = yo(s, e, "ngFalseValue", n.ngFalseValue, !1);
          t.on("change", (function(e) {
            r.$setViewValue(t[0].checked, e && e.type)
          })), r.$render = function() {
            t[0].checked = r.$viewValue
          }, r.$isEmpty = function(e) {
            return !1 === e
          }, r.$formatters.push((function(e) {
            return ue(e, u)
          })), r.$parsers.push((function(e) {
            return e ? u : l
          }))
        },
        hidden: q,
        button: q,
        submit: q,
        reset: q,
        file: q
      };

      function lo(e) {
        e.$formatters.push((function(t) {
          return e.$isEmpty(t) ? t : t.toString()
        }))
      }

      function co(e, t, n, r, i, o) {
        var a, s = p(t[0].type);
        if (!i.android) {
          var u = !1;
          t.on("compositionstart", (function() {
            u = !0
          })), t.on("compositionupdate", (function(e) {
            (j(e.data) || "" === e.data) && (u = !1)
          })), t.on("compositionend", (function() {
            u = !1, l()
          }))
        }
        var l = function(e) {
          if (a && (o.defer.cancel(a), a = null), !u) {
            var i = t.val(),
              l = e && e.type;
            "password" === s || n.ngTrim && "false" === n.ngTrim || (i = ee(i)), (r.$viewValue !== i || "" === i && r.$$hasNativeValidators) && r
              .$setViewValue(i, l)
          }
        };
        if (i.hasEvent("input")) t.on("input", l);
        else {
          var c = function(e, t, n) {
            a || (a = o.defer((function() {
              a = null, t && t.value === n || l(e)
            })))
          };
          t.on("keydown", (function(e) {
            var t = e.keyCode;
            91 === t || 15 < t && t < 19 || 37 <= t && t <= 40 || c(e, this, this.value)
          })), i.hasEvent("paste") && t.on("paste cut drop", c)
        }
        t.on("change", l), so[s] && r.$$hasNativeValidators && s === n.type && t.on("keydown wheel mousedown", (function(e) {
          if (!a) {
            var t = this.validity,
              n = t.badInput,
              r = t.typeMismatch;
            a = o.defer((function() {
              a = null, t.badInput === n && t.typeMismatch === r || l(e)
            }))
          }
        })), r.$render = function() {
          var e = r.$isEmpty(r.$viewValue) ? "" : r.$viewValue;
          t.val() !== e && t.val(e)
        }
      }

      function po(e, t) {
        return function(n, r) {
          var i, o;
          if (B(n)) return n;
          if (H(n)) {
            if ('"' === n.charAt(0) && '"' === n.charAt(n.length - 1) && (n = n.substring(1, n.length - 1)), Ji.test(n)) return new Date(n);
            if (e.lastIndex = 0, i = e.exec(n)) {
              i.shift(), o = r ? {
                yyyy: r.getFullYear(),
                MM: r.getMonth() + 1,
                dd: r.getDate(),
                HH: r.getHours(),
                mm: r.getMinutes(),
                ss: r.getSeconds(),
                sss: r.getMilliseconds() / 1e3
              } : {
                yyyy: 1970,
                MM: 1,
                dd: 1,
                HH: 0,
                mm: 0,
                ss: 0,
                sss: 0
              }, C(i, (function(e, n) {
                n < t.length && (o[t[n]] = +e)
              }));
              var a = new Date(o.yyyy, o.MM - 1, o.dd, o.HH, o.mm, o.ss || 0, 1e3 * o.sss || 0);
              return o.yyyy < 100 && a.setFullYear(o.yyyy), a
            }
          }
          return NaN
        }
      }

      function fo(e, t, n, r) {
        return function(i, o, a, s, u, l, c, d) {
          ho(i, o, a, s, e), co(0, o, a, s, u, l);
          var p, f, h = "time" === e || "datetimelocal" === e;
          if (s.$parsers.push((function(n) {
              return s.$isEmpty(n) ? null : t.test(n) ? w(n, p) : void(s.$$parserName = e)
            })), s.$formatters.push((function(e) {
              if (e && !B(e)) throw Xo("datefmt", "Expected `{0}` to be a date", e);
              if (b(e)) {
                p = e;
                var t = s.$options.getOption("timezone");
                return t && (f = t, p = ye(p, t, !0)),
                  function(e, t) {
                    var n = r;
                    h && H(s.$options.getOption("timeSecondsFormat")) && (n = r.replace("ss.sss", s.$options.getOption("timeSecondsFormat"))
                      .replace(/:$/, ""));
                    var i = c("date")(e, n, t);
                    h && s.$options.getOption("timeStripZeroSeconds") && (i = i.replace(/(?::00)?(?:\.000)?$/, ""));
                    return i
                  }(e, t)
              }
              return p = null, f = null, ""
            })), _(a.min) || a.ngMin) {
            var g = a.min || d(a.ngMin)(i),
              m = y(g);
            s.$validators.min = function(e) {
              return !b(e) || j(m) || n(e) >= m
            }, a.$observe("min", (function(e) {
              e !== g && (m = y(e), g = e, s.$validate())
            }))
          }
          if (_(a.max) || a.ngMax) {
            var v = a.max || d(a.ngMax)(i),
              $ = y(v);
            s.$validators.max = function(e) {
              return !b(e) || j($) || n(e) <= $
            }, a.$observe("max", (function(e) {
              e !== v && ($ = y(e), v = e, s.$validate())
            }))
          }

          function b(e) {
            return e && !(e.getTime && e.getTime() != e.getTime())
          }

          function y(e) {
            return _(e) && !B(e) ? w(e) || void 0 : e
          }

          function w(e, t) {
            var r = s.$options.getOption("timezone");
            f && f !== r && (t = be(t, $e(f)));
            var i = n(e, t);
            return !isNaN(i) && r && (i = ye(i, r)), i
          }
        }
      }

      function ho(e, t, n, r, i) {
        var o = t[0];
        (r.$$hasNativeValidators = V(o.validity)) && r.$parsers.push((function(e) {
          var n = t.prop(c) || {};
          if (!n.badInput && !n.typeMismatch) return e;
          r.$$parserName = i
        }))
      }

      function go(e) {
        e.$parsers.push((function(t) {
          return e.$isEmpty(t) ? null : to.test(t) ? parseFloat(t) : void(e.$$parserName = "number")
        })), e.$formatters.push((function(t) {
          if (!e.$isEmpty(t)) {
            if (!F(t)) throw Xo("numfmt", "Expected `{0}` to be a number", t);
            t = t.toString()
          }
          return t
        }))
      }

      function mo(e) {
        return _(e) && !F(e) && (e = parseFloat(e)), N(e) ? void 0 : e
      }

      function vo(e) {
        return (0 | e) === e
      }

      function $o(e) {
        var t = e.toString(),
          n = t.indexOf(".");
        if (-1 === n) {
          if (-1 < e && e < 1) {
            var r = /e-(\d+)$/.exec(t);
            if (r) return Number(r[1])
          }
          return 0
        }
        return t.length - n - 1
      }

      function bo(e, t, n) {
        var r = Number(e),
          i = !vo(r),
          o = !vo(t),
          a = !vo(n);
        if (i || o || a) {
          var s = i ? $o(r) : 0,
            u = o ? $o(t) : 0,
            l = a ? $o(n) : 0,
            c = Math.max(s, u, l),
            d = Math.pow(10, c);
          r *= d, t *= d, n *= d, i && (r = Math.round(r)), o && (t = Math.round(t)), a && (n = Math.round(n))
        }
        return (r - t) % n == 0
      }

      function yo(e, t, n, r, i) {
        var o;
        if (_(r)) {
          if (!(o = e(r)).constant) throw Xo("constexpr", "Expected constant expression for `{0}`, but saw `{1}`.", n, r);
          return o(t)
        }
        return i
      }
      var wo = ["$browser", "$sniffer", "$filter", "$parse", function(e, t, n, r) {
          return {
            restrict: "E",
            require: ["?ngModel"],
            link: {
              pre: function(i, o, a, s) {
                s[0] && (uo[p(a.type)] || uo.text)(i, o, a, s[0], t, e, n, r)
              }
            }
          }
        }],
        xo = function() {
          var e = {
            configurable: !0,
            enumerable: !1,
            get: function() {
              return this.getAttribute("value") || ""
            },
            set: function(e) {
              this.setAttribute("value", e)
            }
          };
          return {
            restrict: "E",
            priority: 200,
            compile: function(t, n) {
              if ("hidden" === p(n.type)) return {
                pre: function(t, n, r, i) {
                  var o = n[0];
                  o.parentNode && o.parentNode.insertBefore(o, o.nextSibling), Object.defineProperty && Object.defineProperty(o, "value", e)
                }
              }
            }
          }
        },
        Co = /^(true|false|\d+)$/,
        So = function() {
          function e(e, t, n) {
            var r = _(n) ? n : 9 === o ? "" : null;
            e.prop("value", r), t.$set("value", n)
          }
          return {
            restrict: "A",
            priority: 100,
            compile: function(t, n) {
              return Co.test(n.ngValue) ? function(t, n, r) {
                e(n, r, t.$eval(r.ngValue))
              } : function(t, n, r) {
                t.$watch(r.ngValue, (function(t) {
                  e(n, r, t)
                }))
              }
            }
          }
        },
        Ao = ["$compile", function(e) {
          return {
            restrict: "AC",
            compile: function(t) {
              return e.$$addBindingClass(t),
                function(t, n, r) {
                  e.$$addBindingInfo(n, r.ngBind), n = n[0], t.$watch(r.ngBind, (function(e) {
                    n.textContent = Ue(e)
                  }))
                }
            }
          }
        }],
        To = ["$interpolate", "$compile", function(e, t) {
          return {
            compile: function(n) {
              return t.$$addBindingClass(n),
                function(n, r, i) {
                  var o = e(r.attr(i.$attr.ngBindTemplate));
                  t.$$addBindingInfo(r, o.expressions), r = r[0], i.$observe("ngBindTemplate", (function(e) {
                    r.textContent = j(e) ? "" : e
                  }))
                }
            }
          }
        }],
        ko = ["$sce", "$parse", "$compile", function(e, t, n) {
          return {
            restrict: "A",
            compile: function(r, i) {
              var o = t(i.ngBindHtml),
                a = t(i.ngBindHtml, (function(t) {
                  return e.valueOf(t)
                }));
              return n.$$addBindingClass(r),
                function(t, r, i) {
                  n.$$addBindingInfo(r, i.ngBindHtml), t.$watch(a, (function() {
                    var n = o(t);
                    r.html(e.getTrustedHtml(n) || "")
                  }))
                }
            }
          }
        }],
        Do = L({
          restrict: "A",
          require: "ngModel",
          link: function(e, t, n, r) {
            r.$viewChangeListeners.push((function() {
              e.$eval(n.ngChange)
            }))
          }
        });

      function Mo(e, t) {
        var n;
        return e = "ngClass" + e, ["$parse", function(a) {
          return {
            restrict: "AC",
            link: function(s, u, l) {
              var c, d = u.data("$classCounts"),
                p = !0;

              function f(e, t) {
                var n = [];
                return C(e, (function(e) {
                  (t > 0 || d[e]) && (d[e] = (d[e] || 0) + t, d[e] === +(t > 0) && n.push(e))
                })), n.join(" ")
              }
              d || (d = Ve(), u.data("$classCounts", d)), "ngClass" !== e && (n || (n = a("$index", (function(e) {
                return 1 & e
              }))), s.$watch(n, (function(e) {
                e === t ? (n = f(i(n = c), 1), l.$addClass(n)) : function(e) {
                  e = f(i(e), -1), l.$removeClass(e)
                }(c);
                var n;
                p = e
              }))), s.$watch(a(l[e], o), (function(e) {
                p === t && function(e, t) {
                  var n = i(e),
                    o = i(t),
                    a = r(n, o),
                    s = r(o, n),
                    u = f(a, -1),
                    c = f(s, 1);
                  l.$addClass(c), l.$removeClass(u)
                }(c, e);
                c = e
              }))
            }
          }
        }];

        function r(e, t) {
          if (!e || !e.length) return [];
          if (!t || !t.length) return e;
          var n = [];
          e: for (var r = 0; r < e.length; r++) {
            for (var i = e[r], o = 0; o < t.length; o++)
              if (i === t[o]) continue e;
            n.push(i)
          }
          return n
        }

        function i(e) {
          return e && e.split(" ")
        }

        function o(e) {
          if (!e) return e;
          var t = e;
          return G(e) ? t = e.map(o).join(" ") : V(e) ? t = Object.keys(e).filter((function(t) {
            return e[t]
          })).join(" ") : H(e) || (t = e + ""), t
        }
      }
      var Eo = Mo("", !0),
        Oo = Mo("Odd", 0),
        No = Mo("Even", 1),
        Po = _i({
          compile: function(e, t) {
            t.$set("ngCloak", void 0), e.removeClass("ng-cloak")
          }
        }),
        qo = [function() {
          return {
            restrict: "A",
            scope: !0,
            controller: "@",
            priority: 500
          }
        }],
        Io = {},
        Lo = {
          blur: !0,
          focus: !0
        };

      function Ro(e, t, n, r, i, o) {
        return {
          restrict: "A",
          compile: function(a, s) {
            var u = e(s[r]);
            return function(e, r) {
              r.on(i, (function(r) {
                var i = function() {
                  u(e, {
                    $event: r
                  })
                };
                if (t.$$phase)
                  if (o) e.$evalAsync(i);
                  else try {
                    i()
                  } catch (e) {
                    n(e)
                  } else e.$apply(i)
              }))
            }
          }
        }
      }
      C("click dblclick mousedown mouseup mouseover mouseout mousemove mouseenter mouseleave keydown keyup keypress submit focus blur copy cut paste"
        .split(" "), (function(e) {
          var t = wn("ng-" + e);
          Io[t] = ["$parse", "$rootScope", "$exceptionHandler", function(n, r, i) {
            return Ro(n, r, i, t, e, Lo[e])
          }]
        }));
      var jo = ["$animate", "$compile", function(e, t) {
          return {
            multiElement: !0,
            transclude: "element",
            priority: 600,
            terminal: !0,
            restrict: "A",
            $$tlb: !0,
            link: function(n, r, i, o, a) {
              var s, u, l;
              n.$watch(i.ngIf, (function(n) {
                n ? u || a((function(n, o) {
                  u = o, n[n.length++] = t.$$createComment("end ngIf", i.ngIf), s = {
                    clone: n
                  }, e.enter(n, r.parent(), r)
                })) : (l && (l.remove(), l = null), u && (u.$destroy(), u = null), s && (l = _e(s.clone), e.leave(l).done((function(e) {
                  !1 !== e && (l = null)
                })), s = null))
              }))
            }
          }
        }],
        _o = ["$templateRequest", "$anchorScroll", "$animate", function(e, t, n) {
          return {
            restrict: "ECA",
            priority: 400,
            terminal: !0,
            transclude: "element",
            controller: y.noop,
            compile: function(r, i) {
              var o = i.ngInclude || i.src,
                a = i.onload || "",
                s = i.autoscroll;
              return function(r, i, u, l, c) {
                var d, p, f, h = 0,
                  g = function() {
                    p && (p.remove(), p = null), d && (d.$destroy(), d = null), f && (n.leave(f).done((function(e) {
                      !1 !== e && (p = null)
                    })), p = f, f = null)
                  };
                r.$watch(o, (function(o) {
                  var u = function(e) {
                      !1 === e || !_(s) || s && !r.$eval(s) || t()
                    },
                    p = ++h;
                  o ? (e(o, !0).then((function(e) {
                    if (!r.$$destroyed && p === h) {
                      var t = r.$new();
                      l.template = e;
                      var s = c(t, (function(e) {
                        g(), n.enter(e, null, i).done(u)
                      }));
                      f = s, (d = t).$emit("$includeContentLoaded", o), r.$eval(a)
                    }
                  }), (function() {
                    r.$$destroyed || p === h && (g(), r.$emit("$includeContentError", o))
                  })), r.$emit("$includeContentRequested", o)) : (g(), l.template = null)
                }))
              }
            }
          }
        }],
        Vo = ["$compile", function(t) {
          return {
            restrict: "ECA",
            priority: -400,
            require: "ngInclude",
            link: function(n, r, i, o) {
              if (v.call(r[0]).match(/SVG/)) return r.empty(), void t(dt(o.template, e.document).childNodes)(n, (function(e) {
                r.append(e)
              }), {
                futureParentElement: r
              });
              r.html(o.template), t(r.contents())(n)
            }
          }
        }],
        Uo = _i({
          priority: 450,
          compile: function() {
            return {
              pre: function(e, t, n) {
                e.$eval(n.ngInit)
              }
            }
          }
        }),
        Ho = function() {
          return {
            restrict: "A",
            priority: 100,
            require: "ngModel",
            link: function(e, t, n, r) {
              var i = n.ngList || ", ",
                o = "false" !== n.ngTrim,
                a = o ? ee(i) : i;
              r.$parsers.push((function(e) {
                if (!j(e)) {
                  var t = [];
                  return e && C(e.split(a), (function(e) {
                    e && t.push(o ? ee(e) : e)
                  })), t
                }
              })), r.$formatters.push((function(e) {
                if (G(e)) return e.join(i)
              })), r.$isEmpty = function(e) {
                return !e || !e.length
              }
            }
          }
        },
        Fo = "ng-valid",
        Bo = "ng-invalid",
        Go = "ng-pristine",
        Wo = "ng-dirty",
        zo = "ng-untouched",
        Yo = "ng-touched",
        Ko = "ng-empty",
        Qo = "ng-not-empty",
        Xo = i("ngModel");

      function Jo(e, t, n, r, i, o, a, s, u) {
        var l;
        this.$viewValue = Number.NaN, this.$modelValue = Number.NaN, this.$$rawModelValue = void 0, this.$validators = {}, this.$asyncValidators = {},
          this.$parsers = [], this.$formatters = [], this.$viewChangeListeners = [], this.$untouched = !0, this.$touched = !1, this.$pristine = !0, this
          .$dirty = !1, this.$valid = !0, this.$invalid = !1, this.$error = {}, this.$$success = {}, this.$pending = void 0, this.$name = u(n.name || "",
            !1)(e), this.$$parentForm = Hi, this.$options = Zo, this.$$updateEvents = "", this.$$updateEventHandler = this.$$updateEventHandler.bind(
          this), this.$$parsedNgModel = i(n.ngModel), this.$$parsedNgModelAssign = this.$$parsedNgModel.assign, this.$$ngModelGet = this.$$parsedNgModel,
          this.$$ngModelSet = this.$$parsedNgModelAssign, this.$$pendingDebounce = null, this.$$parserValid = void 0, this.$$parserName = "parse", this
          .$$currentValidationRunId = 0, this.$$scope = e, this.$$rootScope = e.$root, this.$$attr = n, this.$$element = r, this.$$animate = o, this
          .$$timeout = a, this.$$parse = i, this.$$q = s, this.$$exceptionHandler = t, Ki(this), (l = this).$$scope.$watch((function(e) {
            var t = l.$$ngModelGet(e);
            return t === l.$modelValue || l.$modelValue != l.$modelValue && t != t || l.$$setModelValue(t), t
          }))
      }
      Jo.$inject = ["$scope", "$exceptionHandler", "$attrs", "$element", "$parse", "$animate", "$timeout", "$q", "$interpolate"], Jo.prototype = {
        $$initGetterSetters: function() {
          if (this.$options.getOption("getterSetter")) {
            var e = this.$$parse(this.$$attr.ngModel + "()"),
              t = this.$$parse(this.$$attr.ngModel + "($$$p)");
            this.$$ngModelGet = function(t) {
              var n = this.$$parsedNgModel(t);
              return z(n) && (n = e(t)), n
            }, this.$$ngModelSet = function(e, n) {
              z(this.$$parsedNgModel(e)) ? t(e, {
                $$$p: n
              }) : this.$$parsedNgModelAssign(e, n)
            }
          } else if (!this.$$parsedNgModel.assign) throw Xo("nonassign", "Expression '{0}' is non-assignable. Element: {1}", this.$$attr.ngModel, we(
            this.$$element))
        },
        $render: q,
        $isEmpty: function(e) {
          return j(e) || "" === e || null === e || e != e
        },
        $$updateEmptyClasses: function(e) {
          this.$isEmpty(e) ? (this.$$animate.removeClass(this.$$element, Qo), this.$$animate.addClass(this.$$element, Ko)) : (this.$$animate
            .removeClass(this.$$element, Ko), this.$$animate.addClass(this.$$element, Qo))
        },
        $setPristine: function() {
          this.$dirty = !1, this.$pristine = !0, this.$$animate.removeClass(this.$$element, Wo), this.$$animate.addClass(this.$$element, Go)
        },
        $setDirty: function() {
          this.$dirty = !0, this.$pristine = !1, this.$$animate.removeClass(this.$$element, Go), this.$$animate.addClass(this.$$element, Wo), this
            .$$parentForm.$setDirty()
        },
        $setUntouched: function() {
          this.$touched = !1, this.$untouched = !0, this.$$animate.setClass(this.$$element, zo, Yo)
        },
        $setTouched: function() {
          this.$touched = !0, this.$untouched = !1, this.$$animate.setClass(this.$$element, Yo, zo)
        },
        $rollbackViewValue: function() {
          this.$$timeout.cancel(this.$$pendingDebounce), this.$viewValue = this.$$lastCommittedViewValue, this.$render()
        },
        $validate: function() {
          if (!N(this.$modelValue)) {
            var e = this.$$lastCommittedViewValue,
              t = this.$$rawModelValue,
              n = this.$valid,
              r = this.$modelValue,
              i = this.$options.getOption("allowInvalid"),
              o = this;
            this.$$runValidators(t, e, (function(e) {
              i || n === e || (o.$modelValue = e ? t : void 0, o.$modelValue !== r && o.$$writeModelToScope())
            }))
          }
        },
        $$runValidators: function(e, t, n) {
          this.$$currentValidationRunId++;
          var r, i, o = this.$$currentValidationRunId,
            a = this;
          (function() {
            var e = a.$$parserName;
            if (!j(a.$$parserValid)) return a.$$parserValid || (C(a.$validators, (function(e, t) {
              s(t, null)
            })), C(a.$asyncValidators, (function(e, t) {
              s(t, null)
            }))), s(e, a.$$parserValid), a.$$parserValid;
            s(e, null);
            return !0
          })() ? ! function() {
            var n = !0;
            if (C(a.$validators, (function(r, i) {
                var o = Boolean(r(e, t));
                n = n && o, s(i, o)
              })), !n) return C(a.$asyncValidators, (function(e, t) {
              s(t, null)
            })), !1;
            return !0
          }() ? u(!1): (r = [], i = !0, C(a.$asyncValidators, (function(n, o) {
            var a = n(e, t);
            if (!J(a)) throw Xo("nopromise", "Expected asynchronous validator to return a promise but got '{0}' instead.", a);
            s(o, void 0), r.push(a.then((function() {
              s(o, !0)
            }), (function() {
              i = !1, s(o, !1)
            })))
          })), r.length ? a.$$q.all(r).then((function() {
            u(i)
          }), q) : u(!0)): u(!1);

          function s(e, t) {
            o === a.$$currentValidationRunId && a.$setValidity(e, t)
          }

          function u(e) {
            o === a.$$currentValidationRunId && n(e)
          }
        },
        $commitViewValue: function() {
          var e = this.$viewValue;
          this.$$timeout.cancel(this.$$pendingDebounce), (this.$$lastCommittedViewValue !== e || "" === e && this.$$hasNativeValidators) && (this
            .$$updateEmptyClasses(e), this.$$lastCommittedViewValue = e, this.$pristine && this.$setDirty(), this.$$parseAndValidate())
        },
        $$parseAndValidate: function() {
          var e = this.$$lastCommittedViewValue,
            t = this;
          if (this.$$parserValid = !j(e) || void 0, this.$setValidity(this.$$parserName, null), this.$$parserName = "parse", this.$$parserValid)
            for (var n = 0; n < this.$parsers.length; n++)
              if (j(e = this.$parsers[n](e))) {
                this.$$parserValid = !1;
                break
              } N(this.$modelValue) && (this.$modelValue = this.$$ngModelGet(this.$$scope));
          var r = this.$modelValue,
            i = this.$options.getOption("allowInvalid");

          function o() {
            t.$modelValue !== r && t.$$writeModelToScope()
          }
          this.$$rawModelValue = e, i && (this.$modelValue = e, o()), this.$$runValidators(e, this.$$lastCommittedViewValue, (function(n) {
            i || (t.$modelValue = n ? e : void 0, o())
          }))
        },
        $$writeModelToScope: function() {
          this.$$ngModelSet(this.$$scope, this.$modelValue), C(this.$viewChangeListeners, (function(e) {
            try {
              e()
            } catch (e) {
              this.$$exceptionHandler(e)
            }
          }), this)
        },
        $setViewValue: function(e, t) {
          this.$viewValue = e, this.$options.getOption("updateOnDefault") && this.$$debounceViewValueCommit(t)
        },
        $$debounceViewValueCommit: function(e) {
          var t = this.$options.getOption("debounce");
          F(t[e]) ? t = t[e] : F(t.default) && -1 === this.$options.getOption("updateOn").indexOf(e) ? t = t.default : F(t["*"]) && (t = t["*"]), this
            .$$timeout.cancel(this.$$pendingDebounce);
          var n = this;
          t > 0 ? this.$$pendingDebounce = this.$$timeout((function() {
            n.$commitViewValue()
          }), t) : this.$$rootScope.$$phase ? this.$commitViewValue() : this.$$scope.$apply((function() {
            n.$commitViewValue()
          }))
        },
        $overrideModelOptions: function(e) {
          this.$options = this.$options.createChild(e), this.$$setUpdateOnEvents()
        },
        $processModelValue: function() {
          var e = this.$$format();
          this.$viewValue !== e && (this.$$updateEmptyClasses(e), this.$viewValue = this.$$lastCommittedViewValue = e, this.$render(), this
            .$$runValidators(this.$modelValue, this.$viewValue, q))
        },
        $$format: function() {
          for (var e = this.$formatters, t = e.length, n = this.$modelValue; t--;) n = e[t](n);
          return n
        },
        $$setModelValue: function(e) {
          this.$modelValue = this.$$rawModelValue = e, this.$$parserValid = void 0, this.$processModelValue()
        },
        $$setUpdateOnEvents: function() {
          this.$$updateEvents && this.$$element.off(this.$$updateEvents, this.$$updateEventHandler), this.$$updateEvents = this.$options.getOption(
            "updateOn"), this.$$updateEvents && this.$$element.on(this.$$updateEvents, this.$$updateEventHandler)
        },
        $$updateEventHandler: function(e) {
          this.$$debounceViewValueCommit(e && e.type)
        }
      }, Qi({
        clazz: Jo,
        set: function(e, t) {
          e[t] = !0
        },
        unset: function(e, t) {
          delete e[t]
        }
      });
      var Zo, ea = ["$rootScope", function(e) {
          return {
            restrict: "A",
            require: ["ngModel", "^?form", "^?ngModelOptions"],
            controller: Jo,
            priority: 1,
            compile: function(t) {
              return t.addClass(Go).addClass(zo).addClass(Fo), {
                pre: function(e, t, n, r) {
                  var i = r[0],
                    o = r[1] || i.$$parentForm,
                    a = r[2];
                  a && (i.$options = a.$options), i.$$initGetterSetters(), o.$addControl(i), n.$observe("name", (function(e) {
                    i.$name !== e && i.$$parentForm.$$renameControl(i, e)
                  })), e.$on("$destroy", (function() {
                    i.$$parentForm.$removeControl(i)
                  }))
                },
                post: function(t, n, r, i) {
                  var o = i[0];

                  function a() {
                    o.$setTouched()
                  }
                  o.$$setUpdateOnEvents(), n.on("blur", (function() {
                    o.$touched || (e.$$phase ? t.$evalAsync(a) : t.$apply(a))
                  }))
                }
              }
            }
          }
        }],
        ta = /(\s+|^)default(\s+|$)/;

      function na(e) {
        this.$$options = e
      }
      na.prototype = {
        getOption: function(e) {
          return this.$$options[e]
        },
        createChild: function(e) {
          var t = !1;
          return C(e = M({}, e), (function(n, r) {
            "$inherit" === n ? "*" === r ? t = !0 : (e[r] = this.$$options[r], "updateOn" === r && (e.updateOnDefault = this.$$options
              .updateOnDefault)) : "updateOn" === r && (e.updateOnDefault = !1, e[r] = ee(n.replace(ta, (function() {
              return e.updateOnDefault = !0, " "
            }))))
          }), this), t && (delete e["*"], ia(e, this.$$options)), ia(e, Zo.$$options), new na(e)
        }
      }, Zo = new na({
        updateOn: "",
        updateOnDefault: !0,
        debounce: 0,
        getterSetter: !1,
        allowInvalid: !1,
        timezone: null
      });
      var ra = function() {
        function e(e, t) {
          this.$$attrs = e, this.$$scope = t
        }
        return e.$inject = ["$attrs", "$scope"], e.prototype = {
          $onInit: function() {
            var e = this.parentCtrl ? this.parentCtrl.$options : Zo,
              t = this.$$scope.$eval(this.$$attrs.ngModelOptions);
            this.$options = e.createChild(t)
          }
        }, {
          restrict: "A",
          priority: 10,
          require: {
            parentCtrl: "?^^ngModelOptions"
          },
          bindToController: !0,
          controller: e
        }
      };

      function ia(e, t) {
        C(t, (function(t, n) {
          _(e[n]) || (e[n] = t)
        }))
      }
      var oa = _i({
          terminal: !0,
          priority: 1e3
        }),
        aa = i("ngOptions"),
        sa =
        /^\s*([\s\S]+?)(?:\s+as\s+([\s\S]+?))?(?:\s+group\s+by\s+([\s\S]+?))?(?:\s+disable\s+when\s+([\s\S]+?))?\s+for\s+(?:([$\w][$\w]*)|(?:\(\s*([$\w][$\w]*)\s*,\s*([$\w][$\w]*)\s*\)))\s+in\s+([\s\S]+?)(?:\s+track\s+by\s+([\s\S]+?))?$/,
        ua = ["$compile", "$document", "$parse", function(t, n, r) {
          var i = e.document.createElement("option"),
            o = e.document.createElement("optgroup");
          return {
            restrict: "A",
            terminal: !0,
            require: ["select", "ngModel"],
            link: {
              pre: function(e, t, n, r) {
                r[0].registerOption = q
              },
              post: function(e, s, u, l) {
                for (var c = l[0], d = l[1], p = u.multiple, f = 0, h = s.children(), g = h.length; f < g; f++)
                  if ("" === h[f].value) {
                    c.hasEmptyOption = !0, c.emptyOption = h.eq(f);
                    break
                  } s.empty();
                var m, v = !!c.emptyOption;
                a(i.cloneNode(!1)).val("?");
                var $ = function(e, t, n) {
                    var i = e.match(sa);
                    if (!i) throw aa("iexp",
                      "Expected expression in form of '_select_ (as _label_)? for (_key_,)?_value_ in _collection_' but got '{0}'. Element: {1}",
                      e, we(t));
                    var o = i[5] || i[7],
                      a = i[6],
                      s = / as /.test(i[0]) && i[1],
                      u = i[9],
                      l = r(i[2] ? i[1] : o),
                      c = s && r(s) || l,
                      d = u && r(u),
                      p = u ? function(e, t) {
                        return d(n, t)
                      } : function(e) {
                        return _t(e)
                      },
                      f = function(e, t) {
                        return p(e, b(e, t))
                      },
                      h = r(i[2] || i[1]),
                      g = r(i[3] || ""),
                      m = r(i[4] || ""),
                      v = r(i[8]),
                      $ = {},
                      b = a ? function(e, t) {
                        return $[a] = t, $[o] = e, $
                      } : function(e) {
                        return $[o] = e, $
                      };

                    function y(e, t, n, r, i) {
                      this.selectValue = e, this.viewValue = t, this.label = n, this.group = r, this.disabled = i
                    }

                    function w(e) {
                      var t;
                      if (!a && x(e)) t = e;
                      else
                        for (var n in t = [], e) e.hasOwnProperty(n) && "$" !== n.charAt(0) && t.push(n);
                      return t
                    }
                    return {
                      trackBy: u,
                      getTrackByValue: f,
                      getWatchables: r(v, (function(e) {
                        for (var t = [], r = w(e = e || []), o = r.length, a = 0; a < o; a++) {
                          var s = e === r ? a : r[a],
                            u = e[s],
                            l = b(u, s),
                            c = p(u, l);
                          if (t.push(c), i[2] || i[1]) {
                            var d = h(n, l);
                            t.push(d)
                          }
                          if (i[4]) {
                            var f = m(n, l);
                            t.push(f)
                          }
                        }
                        return t
                      })),
                      getOptions: function() {
                        for (var e = [], t = {}, r = v(n) || [], i = w(r), o = i.length, a = 0; a < o; a++) {
                          var s = r === i ? a : i[a],
                            l = r[s],
                            d = b(l, s),
                            $ = c(n, d),
                            x = p($, d),
                            C = new y(x, $, h(n, d), g(n, d), m(n, d));
                          e.push(C), t[x] = C
                        }
                        return {
                          items: e,
                          selectValueMap: t,
                          getOptionFromViewValue: function(e) {
                            return t[f(e)]
                          },
                          getViewValueFromOption: function(e) {
                            return u ? ae(e.viewValue) : e.viewValue
                          }
                        }
                      }
                    }
                  }(u.ngOptions, s, e),
                  b = n[0].createDocumentFragment();

                function y(e, t) {
                  var n = i.cloneNode(!1);
                  t.appendChild(n),
                    function(e, t) {
                      e.element = t, t.disabled = e.disabled, e.label !== t.label && (t.label = e.label, t.textContent = e.label);
                      t.value = e.selectValue
                    }(e, n)
                }

                function w(e) {
                  var t = m.getOptionFromViewValue(e),
                    n = t && t.element;
                  return n && !n.selected && (n.selected = !0), t
                }
                c.generateUnknownOptionValue = function(e) {
                  return "?"
                }, p ? (c.writeValue = function(e) {
                  if (m) {
                    var t = e && e.map(w) || [];
                    m.items.forEach((function(e) {
                      e.element.selected && !ie(t, e) && (e.element.selected = !1)
                    }))
                  }
                }, c.readValue = function() {
                  var e = s.val() || [],
                    t = [];
                  return C(e, (function(e) {
                    var n = m.selectValueMap[e];
                    n && !n.disabled && t.push(m.getViewValueFromOption(n))
                  })), t
                }, $.trackBy && e.$watchCollection((function() {
                  if (G(d.$viewValue)) return d.$viewValue.map((function(e) {
                    return $.getTrackByValue(e)
                  }))
                }), (function() {
                  d.$render()
                }))) : (c.writeValue = function(e) {
                  if (m) {
                    var t = s[0].options[s[0].selectedIndex],
                      n = m.getOptionFromViewValue(e);
                    t && t.removeAttribute("selected"), n ? (s[0].value !== n.selectValue && (c.removeUnknownOption(), s[0].value = n.selectValue,
                      n.element.selected = !0), n.element.setAttribute("selected", "selected")) : c.selectUnknownOrEmptyOption(e)
                  }
                }, c.readValue = function() {
                  var e = m.selectValueMap[s.val()];
                  return e && !e.disabled ? (c.unselectEmptyOption(), c.removeUnknownOption(), m.getViewValueFromOption(e)) : null
                }, $.trackBy && e.$watch((function() {
                  return $.getTrackByValue(d.$viewValue)
                }), (function() {
                  d.$render()
                }))), v && (t(c.emptyOption)(e), s.prepend(c.emptyOption), 8 === c.emptyOption[0].nodeType ? (c.hasEmptyOption = !1, c
                  .registerOption = function(e, t) {
                    "" === t.val() && (c.hasEmptyOption = !0, c.emptyOption = t, c.emptyOption.removeClass("ng-scope"), d.$render(), t.on(
                      "$destroy", (function() {
                        var e = c.$isEmptyOptionSelected();
                        c.hasEmptyOption = !1, c.emptyOption = void 0, e && d.$render()
                      })))
                  }) : c.emptyOption.removeClass("ng-scope")), e.$watchCollection($.getWatchables, (function() {
                  var e = m && c.readValue();
                  if (m)
                    for (var t = m.items.length - 1; t >= 0; t--) {
                      var n = m.items[t];
                      _(n.group) ? Mt(n.element.parentNode) : Mt(n.element)
                    }
                  m = $.getOptions();
                  var r = {};
                  if (m.items.forEach((function(e) {
                      var t;
                      _(e.group) ? ((t = r[e.group]) || (t = o.cloneNode(!1), b.appendChild(t), t.label = null === e.group ? "null" : e
                        .group, r[e.group] = t), y(e, t)) : y(e, b)
                    })), s[0].appendChild(b), d.$render(), !d.$isEmpty(e)) {
                    var i = c.readValue();
                    ($.trackBy || p ? ue(e, i) : e === i) || (d.$setViewValue(i), d.$render())
                  }
                }))
              }
            }
          }
        }],
        la = ["$locale", "$interpolate", "$log", function(e, t, n) {
          var r = /{}/g,
            i = /^when(Minus)?(.+)$/;
          return {
            link: function(o, a, s) {
              var u, l = s.count,
                c = s.$attr.when && a.attr(s.$attr.when),
                d = s.offset || 0,
                f = o.$eval(c) || {},
                h = {},
                g = t.startSymbol(),
                m = t.endSymbol(),
                v = g + l + "-" + d + m,
                $ = y.noop;

              function b(e) {
                a.text(e || "")
              }
              C(s, (function(e, t) {
                var n = i.exec(t);
                if (n) {
                  var r = (n[1] ? "-" : "") + p(n[2]);
                  f[r] = a.attr(s.$attr[t])
                }
              })), C(f, (function(e, n) {
                h[n] = t(e.replace(r, v))
              })), o.$watch(l, (function(t) {
                var r = parseFloat(t),
                  i = N(r);
                if (i || r in f || (r = e.pluralCat(r - d)), !(r === u || i && N(u))) {
                  $();
                  var a = h[r];
                  j(a) ? (null != t && n.debug("ngPluralize: no rule defined for '" + r + "' in " + c), $ = q, b()) : $ = o.$watch(a, b), u = r
                }
              }))
            }
          }
        }],
        ca = i("ngRef"),
        da = ["$parse", function(e) {
          return {
            priority: -1,
            restrict: "A",
            compile: function(t, n) {
              var r = wn(re(t)),
                i = e(n.ngRef),
                o = i.assign || function() {
                  throw ca("nonassign", 'Expression in ngRef="{0}" is non-assignable!', n.ngRef)
                };
              return function(e, t, a) {
                var s;
                if (a.hasOwnProperty("ngRefRead")) {
                  if ("$element" === a.ngRefRead) s = t;
                  else if (!(s = t.data("$" + a.ngRefRead + "Controller"))) throw ca("noctrl",
                    'The controller for ngRefRead="{0}" could not be found on ngRef="{1}"', a.ngRefRead, n.ngRef)
                } else s = t.data("$" + r + "Controller");
                o(e, s = s || t), t.on("$destroy", (function() {
                  i(e) === s && o(e, null)
                }))
              }
            }
          }
        }],
        pa = ["$parse", "$animate", "$compile", function(e, t, n) {
          var r = "$$NG_REMOVED",
            o = i("ngRepeat"),
            a = function(e, t, n, r, i, o, a) {
              e[n] = r, i && (e[i] = o), e.$index = t, e.$first = 0 === t, e.$last = t === a - 1, e.$middle = !(e.$first || e.$last), e.$odd = !(e
                .$even = 0 == (1 & t))
            },
            s = function(e) {
              return e.clone[0]
            },
            u = function(e) {
              return e.clone[e.clone.length - 1]
            },
            l = function(e, t, n) {
              return _t(n)
            },
            c = function(e, t) {
              return t
            };
          return {
            restrict: "A",
            multiElement: !0,
            transclude: "element",
            priority: 1e3,
            terminal: !0,
            $$tlb: !0,
            compile: function(i, p) {
              var f = p.ngRepeat,
                h = n.$$createComment("end ngRepeat", f),
                g = f.match(/^\s*([\s\S]+?)\s+in\s+([\s\S]+?)(?:\s+as\s+([\s\S]+?))?(?:\s+track\s+by\s+([\s\S]+?))?\s*$/);
              if (!g) throw o("iexp", "Expected expression in form of '_item_ in _collection_[ track by _id_]' but got '{0}'.", f);
              var m = g[1],
                v = g[2],
                $ = g[3],
                b = g[4];
              if (!(g = m.match(/^(?:(\s*[$\w]+)|\(\s*([$\w]+)\s*,\s*([$\w]+)\s*\))$/))) throw o("iidexp",
                "'_item_' in '_item_ in _collection_' should be an identifier or '(_key_, _value_)' expression, but got '{0}'.", m);
              var y, w = g[3] || g[1],
                S = g[2];
              if ($ && (!/^[$a-zA-Z_][$a-zA-Z0-9_]*$/.test($) ||
                  /^(null|undefined|this|\$index|\$first|\$middle|\$last|\$even|\$odd|\$parent|\$root|\$id)$/.test($))) throw o("badident",
                "alias '{0}' is invalid --- must be a valid JS identifier which is not a reserved name.", $);
              if (b) {
                var A = {
                    $id: _t
                  },
                  T = e(b);
                y = function(e, t, n, r) {
                  return S && (A[S] = t), A[w] = n, A.$index = r, T(e, A)
                }
              }
              return function(e, n, i, p, g) {
                var m = Ve();
                e.$watchCollection(v, (function(i) {
                  var p, v, b, T, k, D, M, E, O, N, P, q, I = n[0],
                    L = Ve();
                  if ($ && (e[$] = i), x(i)) O = i, E = y || l;
                  else
                    for (var R in E = y || c, O = [], i) d.call(i, R) && "$" !== R.charAt(0) && O.push(R);
                  for (T = O.length, P = new Array(T), p = 0; p < T; p++)
                    if (k = i === O ? p : O[p], D = i[k], M = E(e, k, D, p), m[M]) N = m[M], delete m[M], L[M] = N, P[p] = N;
                    else {
                      if (L[M]) throw C(P, (function(e) {
                        e && e.scope && (m[e.id] = e)
                      })), o("dupes",
                        "Duplicates in a repeater are not allowed. Use 'track by' expression to specify unique keys. Repeater: {0}, Duplicate key: {1}, Duplicate value: {2}",
                        f, M, D);
                      P[p] = {
                        id: M,
                        scope: void 0,
                        clone: void 0
                      }, L[M] = !0
                    } for (var j in A && (A[w] = void 0), m) {
                    if (q = _e((N = m[j]).clone), t.leave(q), q[0].parentNode)
                      for (p = 0, v = q.length; p < v; p++) q[p][r] = !0;
                    N.scope.$destroy()
                  }
                  for (p = 0; p < T; p++)
                    if (k = i === O ? p : O[p], D = i[k], (N = P[p]).scope) {
                      b = I;
                      do {
                        b = b.nextSibling
                      } while (b && b[r]);
                      s(N) !== b && t.move(_e(N.clone), null, I), I = u(N), a(N.scope, p, w, D, S, k, T)
                    } else g((function(e, n) {
                      N.scope = n;
                      var r = h.cloneNode(!1);
                      e[e.length++] = r, t.enter(e, null, I), I = r, N.clone = e, L[N.id] = N, a(N.scope, p, w, D, S, k, T)
                    }));
                  m = L
                }))
              }
            }
          }
        }],
        fa = "ng-hide",
        ha = "ng-hide-animate",
        ga = ["$animate", function(e) {
          return {
            restrict: "A",
            multiElement: !0,
            link: function(t, n, r) {
              t.$watch(r.ngShow, (function(t) {
                e[t ? "removeClass" : "addClass"](n, fa, {
                  tempClasses: ha
                })
              }))
            }
          }
        }],
        ma = ["$animate", function(e) {
          return {
            restrict: "A",
            multiElement: !0,
            link: function(t, n, r) {
              t.$watch(r.ngHide, (function(t) {
                e[t ? "addClass" : "removeClass"](n, fa, {
                  tempClasses: ha
                })
              }))
            }
          }
        }],
        va = _i((function(e, t, n) {
          e.$watchCollection(n.ngStyle, (function(e, n) {
            n && e !== n && C(n, (function(e, n) {
              t.css(n, "")
            })), e && t.css(e)
          }))
        })),
        $a = ["$animate", "$compile", function(e, t) {
          return {
            require: "ngSwitch",
            controller: ["$scope", function() {
              this.cases = {}
            }],
            link: function(n, r, i, o) {
              var a = i.ngSwitch || i.on,
                s = [],
                u = [],
                l = [],
                c = [],
                d = function(e, t) {
                  return function(n) {
                    !1 !== n && e.splice(t, 1)
                  }
                };
              n.$watch(a, (function(n) {
                for (var r, i; l.length;) e.cancel(l.pop());
                for (r = 0, i = c.length; r < i; ++r) {
                  var a = _e(u[r].clone);
                  c[r].$destroy(), (l[r] = e.leave(a)).done(d(l, r))
                }
                u.length = 0, c.length = 0, (s = o.cases["!" + n] || o.cases["?"]) && C(s, (function(n) {
                  n.transclude((function(r, i) {
                    c.push(i);
                    var o = n.element;
                    r[r.length++] = t.$$createComment("end ngSwitchWhen");
                    var a = {
                      clone: r
                    };
                    u.push(a), e.enter(r, o.parent(), o)
                  }))
                }))
              }))
            }
          }
        }],
        ba = _i({
          transclude: "element",
          priority: 1200,
          require: "^ngSwitch",
          multiElement: !0,
          link: function(e, t, n, r, i) {
            C(n.ngSwitchWhen.split(n.ngSwitchWhenSeparator).sort().filter((function(e, t, n) {
              return n[t - 1] !== e
            })), (function(e) {
              r.cases["!" + e] = r.cases["!" + e] || [], r.cases["!" + e].push({
                transclude: i,
                element: t
              })
            }))
          }
        }),
        ya = _i({
          transclude: "element",
          priority: 1200,
          require: "^ngSwitch",
          multiElement: !0,
          link: function(e, t, n, r, i) {
            r.cases["?"] = r.cases["?"] || [], r.cases["?"].push({
              transclude: i,
              element: t
            })
          }
        }),
        wa = i("ngTransclude"),
        xa = ["$compile", function(e) {
          return {
            restrict: "EAC",
            compile: function(t) {
              var n = e(t.contents());
              return t.empty(),
                function(e, t, r, i, o) {
                  if (!o) throw wa("orphan",
                    "Illegal use of ngTransclude directive in the template! No parent directive that requires a transclusion found. Element: {0}",
                    we(t));
                  r.ngTransclude === r.$attr.ngTransclude && (r.ngTransclude = "");
                  var a = r.ngTransclude || r.ngTranscludeSlot;

                  function s() {
                    n(e, (function(e) {
                      t.append(e)
                    }))
                  }
                  o((function(e, n) {
                    e.length && function(e) {
                      for (var t = 0, n = e.length; t < n; t++) {
                        var r = e[t];
                        if (r.nodeType !== He || r.nodeValue.trim()) return !0
                      }
                    }(e) ? t.append(e) : (s(), n.$destroy())
                  }), null, a), a && !o.isSlotFilled(a) && s()
                }
            }
          }
        }],
        Ca = ["$templateCache", function(e) {
          return {
            restrict: "E",
            terminal: !0,
            compile: function(t, n) {
              if ("text/ng-template" === n.type) {
                var r = n.id,
                  i = t[0].text;
                e.put(r, i)
              }
            }
          }
        }],
        Sa = {
          $setViewValue: q,
          $render: q
        };

      function Aa(e, t) {
        e.prop("selected", t), e.attr("selected", t)
      }
      var Ta = ["$element", "$scope", function(t, n) {
          var r = this,
            i = new Ht;
          r.selectValueMap = {}, r.ngModelCtrl = Sa, r.multiple = !1, r.unknownOption = a(e.document.createElement("option")), r.hasEmptyOption = !1, r
            .emptyOption = void 0, r.renderUnknownOption = function(e) {
              var n = r.generateUnknownOptionValue(e);
              r.unknownOption.val(n), t.prepend(r.unknownOption), Aa(r.unknownOption, !0), t.val(n)
            }, r.updateUnknownOption = function(e) {
              var n = r.generateUnknownOptionValue(e);
              r.unknownOption.val(n), Aa(r.unknownOption, !0), t.val(n)
            }, r.generateUnknownOptionValue = function(e) {
              return "? " + _t(e) + " ?"
            }, r.removeUnknownOption = function() {
              r.unknownOption.parent() && r.unknownOption.remove()
            }, r.selectEmptyOption = function() {
              r.emptyOption && (t.val(""), Aa(r.emptyOption, !0))
            }, r.unselectEmptyOption = function() {
              r.hasEmptyOption && Aa(r.emptyOption, !1)
            }, n.$on("$destroy", (function() {
              r.renderUnknownOption = q
            })), r.readValue = function() {
              var e = t.val(),
                n = e in r.selectValueMap ? r.selectValueMap[e] : e;
              return r.hasOption(n) ? n : null
            }, r.writeValue = function(e) {
              var n = t[0].options[t[0].selectedIndex];
              if (n && Aa(a(n), !1), r.hasOption(e)) {
                r.removeUnknownOption();
                var i = _t(e);
                t.val(i in r.selectValueMap ? i : e);
                var o = t[0].options[t[0].selectedIndex];
                Aa(a(o), !0)
              } else r.selectUnknownOrEmptyOption(e)
            }, r.addOption = function(e, t) {
              if (8 !== t[0].nodeType) {
                je(e, '"option value"'), "" === e && (r.hasEmptyOption = !0, r.emptyOption = t);
                var n = i.get(e) || 0;
                i.set(e, n + 1), s()
              }
            }, r.removeOption = function(e) {
              var t = i.get(e);
              t && (1 === t ? (i.delete(e), "" === e && (r.hasEmptyOption = !1, r.emptyOption = void 0)) : i.set(e, t - 1))
            }, r.hasOption = function(e) {
              return !!i.get(e)
            }, r.$hasEmptyOption = function() {
              return r.hasEmptyOption
            }, r.$isUnknownOptionSelected = function() {
              return t[0].options[0] === r.unknownOption[0]
            }, r.$isEmptyOptionSelected = function() {
              return r.hasEmptyOption && t[0].options[t[0].selectedIndex] === r.emptyOption[0]
            }, r.selectUnknownOrEmptyOption = function(e) {
              null == e && r.emptyOption ? (r.removeUnknownOption(), r.selectEmptyOption()) : r.unknownOption.parent().length ? r.updateUnknownOption(
                e) : r.renderUnknownOption(e)
            };
          var o = !1;

          function s() {
            o || (o = !0, n.$$postDigest((function() {
              o = !1, r.ngModelCtrl.$render()
            })))
          }
          var u = !1;

          function l(e) {
            u || (u = !0, n.$$postDigest((function() {
              n.$$destroyed || (u = !1, r.ngModelCtrl.$setViewValue(r.readValue()), e && r.ngModelCtrl.$render())
            })))
          }
          r.registerOption = function(e, t, n, i, o) {
            var a, u;
            n.$attr.ngValue ? n.$observe("value", (function(e) {
              var n, i = t.prop("selected");
              _(u) && (r.removeOption(a), delete r.selectValueMap[u], n = !0), u = _t(e), a = e, r.selectValueMap[u] = e, r.addOption(e, t), t
                .attr("value", u), n && i && l()
            })) : i ? n.$observe("value", (function(e) {
              var n;
              r.readValue();
              var i = t.prop("selected");
              _(a) && (r.removeOption(a), n = !0), a = e, r.addOption(e, t), n && i && l()
            })) : o ? e.$watch(o, (function(e, i) {
              n.$set("value", e);
              var o = t.prop("selected");
              i !== e && r.removeOption(i), r.addOption(e, t), i && o && l()
            })) : r.addOption(n.value, t);
            n.$observe("disabled", (function(e) {
              ("true" === e || e && t.prop("selected")) && (r.multiple ? l(!0) : (r.ngModelCtrl.$setViewValue(null), r.ngModelCtrl.$render()))
            })), t.on("$destroy", (function() {
              var e = r.readValue(),
                t = n.value;
              r.removeOption(t), s(), (r.multiple && e && -1 !== e.indexOf(t) || e === t) && l(!0)
            }))
          }
        }],
        ka = function() {
          return {
            restrict: "E",
            require: ["select", "?ngModel"],
            controller: Ta,
            priority: 1,
            link: {
              pre: function(e, t, n, r) {
                var i = r[0],
                  o = r[1];
                if (!o) return void(i.registerOption = q);
                if (i.ngModelCtrl = o, t.on("change", (function() {
                    i.removeUnknownOption(), e.$apply((function() {
                      o.$setViewValue(i.readValue())
                    }))
                  })), n.multiple) {
                  i.multiple = !0, i.readValue = function() {
                    var e = [];
                    return C(t.find("option"), (function(t) {
                      if (t.selected && !t.disabled) {
                        var n = t.value;
                        e.push(n in i.selectValueMap ? i.selectValueMap[n] : n)
                      }
                    })), e
                  }, i.writeValue = function(e) {
                    C(t.find("option"), (function(t) {
                      var n = !!e && (ie(e, t.value) || ie(e, i.selectValueMap[t.value]));
                      n !== t.selected && Aa(a(t), n)
                    }))
                  };
                  var s, u = NaN;
                  e.$watch((function() {
                    u !== o.$viewValue || ue(s, o.$viewValue) || (s = Fe(o.$viewValue), o.$render()), u = o.$viewValue
                  })), o.$isEmpty = function(e) {
                    return !e || 0 === e.length
                  }
                }
              },
              post: function(e, t, n, r) {
                var i = r[1];
                if (!i) return;
                var o = r[0];
                i.$render = function() {
                  o.writeValue(i.$viewValue)
                }
              }
            }
          }
        },
        Da = ["$interpolate", function(e) {
          return {
            restrict: "E",
            priority: 100,
            compile: function(t, n) {
              var r, i;
              return _(n.ngValue) || (_(n.value) ? r = e(n.value, !0) : (i = e(t.text(), !0)) || n.$set("value", t.text())),
                function(e, t, n) {
                  var o = "$selectController",
                    a = t.parent(),
                    s = a.data(o) || a.parent().data(o);
                  s && s.registerOption(e, t, n, r, i)
                }
            }
          }
        }],
        Ma = ["$parse", function(e) {
          return {
            restrict: "A",
            require: "?ngModel",
            link: function(t, n, r, i) {
              if (i) {
                var o = r.hasOwnProperty("required") || e(r.ngRequired)(t);
                r.ngRequired || (r.required = !0), i.$validators.required = function(e, t) {
                  return !o || !i.$isEmpty(t)
                }, r.$observe("required", (function(e) {
                  o !== e && (o = e, i.$validate())
                }))
              }
            }
          }
        }],
        Ea = ["$parse", function(e) {
          return {
            restrict: "A",
            require: "?ngModel",
            compile: function(t, n) {
              var r, i;
              return n.ngPattern && (r = n.ngPattern, i = "/" === n.ngPattern.charAt(0) && l.test(n.ngPattern) ? function() {
                  return n.ngPattern
                } : e(n.ngPattern)),
                function(e, t, n, o) {
                  if (o) {
                    var a = n.pattern;
                    n.ngPattern ? a = i(e) : r = n.pattern;
                    var s = Pa(a, r, t);
                    n.$observe("pattern", (function(e) {
                      var n = s;
                      s = Pa(e, r, t), (n && n.toString()) !== (s && s.toString()) && o.$validate()
                    })), o.$validators.pattern = function(e, t) {
                      return o.$isEmpty(t) || j(s) || s.test(t)
                    }
                  }
                }
            }
          }
        }],
        Oa = ["$parse", function(e) {
          return {
            restrict: "A",
            require: "?ngModel",
            link: function(t, n, r, i) {
              if (i) {
                var o = r.maxlength || e(r.ngMaxlength)(t),
                  a = qa(o);
                r.$observe("maxlength", (function(e) {
                  o !== e && (a = qa(e), o = e, i.$validate())
                })), i.$validators.maxlength = function(e, t) {
                  return a < 0 || i.$isEmpty(t) || t.length <= a
                }
              }
            }
          }
        }],
        Na = ["$parse", function(e) {
          return {
            restrict: "A",
            require: "?ngModel",
            link: function(t, n, r, i) {
              if (i) {
                var o = r.minlength || e(r.ngMinlength)(t),
                  a = qa(o) || -1;
                r.$observe("minlength", (function(e) {
                  o !== e && (a = qa(e) || -1, o = e, i.$validate())
                })), i.$validators.minlength = function(e, t) {
                  return i.$isEmpty(t) || t.length >= a
                }
              }
            }
          }
        }];

      function Pa(e, t, n) {
        if (e) {
          if (H(e) && (e = new RegExp("^" + e + "$")), !e.test) throw i("ngPattern")("noregexp", "Expected {0} to be a RegExp but was {1}. Element: {2}",
            t, e, we(n));
          return e
        }
      }

      function qa(e) {
        var t = O(e);
        return N(t) ? -1 : t
      }
      e.angular.bootstrap ? e.console && console.log("WARNING: Tried to load AngularJS more than once.") : (! function() {
        var t;
        if (!qe) {
          var n = ce();
          (s = j(n) ? e.jQuery : n ? e[n] : void 0) && s.fn.on ? (a = s, M(s.fn, {
            scope: Ot.scope,
            isolateScope: Ot.isolateScope,
            controller: Ot.controller,
            injector: Ot.injector,
            inheritedData: Ot.inheritedData
          })) : a = ft, t = a.cleanData, a.cleanData = function(e) {
            for (var n, r, i = 0; null != (r = e[i]); i++)(n = (a._data(r) || {}).events) && n.$destroy && a(r).triggerHandler("$destroy");
            t(e)
          }, y.element = a, qe = !0
        }
      }(), function(t) {
        M(t, {
          errorHandlingConfig: n,
          bootstrap: Me,
          copy: ae,
          extend: M,
          merge: E,
          equals: ue,
          element: a,
          forEach: C,
          injector: Jt,
          noop: q,
          bind: fe,
          toJson: ge,
          fromJson: me,
          identity: I,
          isUndefined: j,
          isDefined: _,
          isString: H,
          isFunction: z,
          isObject: V,
          isNumber: F,
          isElement: ne,
          isArray: G,
          version: Ge,
          isDate: B,
          callbacks: {
            $$counter: 0
          },
          getTestability: Oe,
          reloadWithDebugInfo: Ee,
          UNSAFE_restoreLegacyJqLiteXHTMLReplacement: Ie,
          $$minErr: i,
          $$csp: le,
          $$encodeUriSegment: Se,
          $$encodeUriQuery: Ae,
          $$lowercase: p,
          $$stringify: Ue,
          $$uppercase: f
        }), (u = function(e) {
          var t = i("$injector"),
            n = i("ng");

          function r(e, t, n) {
            return e[t] || (e[t] = n())
          }
          var o = r(e, "angular", Object);
          return o.$$minErr = o.$$minErr || i, r(o, "module", (function() {
            var e = {};
            return function(i, o, a) {
              var s = {};
              return function(e, t) {
                if ("hasOwnProperty" === e) throw n("badname", "hasOwnProperty is not a valid {0} name", "module")
              }(i), o && e.hasOwnProperty(i) && (e[i] = null), r(e, i, (function() {
                if (!o) throw t("nomod",
                  "Module '{0}' is not available! You either misspelled the module name or forgot to load it. If registering a module ensure that you specify the dependencies as the second argument.",
                  i);
                var e = [],
                  r = [],
                  u = [],
                  l = d("$injector", "invoke", "push", r),
                  c = {
                    _invokeQueue: e,
                    _configBlocks: r,
                    _runBlocks: u,
                    info: function(e) {
                      if (_(e)) {
                        if (!V(e)) throw n("aobj", "Argument '{0}' must be an object", "value");
                        return s = e, this
                      }
                      return s
                    },
                    requires: o,
                    name: i,
                    provider: p("$provide", "provider"),
                    factory: p("$provide", "factory"),
                    service: p("$provide", "service"),
                    value: d("$provide", "value"),
                    constant: d("$provide", "constant", "unshift"),
                    decorator: p("$provide", "decorator", r),
                    animation: p("$animateProvider", "register"),
                    filter: p("$filterProvider", "register"),
                    controller: p("$controllerProvider", "register"),
                    directive: p("$compileProvider", "directive"),
                    component: p("$compileProvider", "component"),
                    config: l,
                    run: function(e) {
                      return u.push(e), this
                    }
                  };
                return a && l(a), c;

                function d(t, n, r, i) {
                  return i || (i = e),
                    function() {
                      return i[r || "push"]([t, n, arguments]), c
                    }
                }

                function p(t, n, r) {
                  return r || (r = e),
                    function(e, o) {
                      return o && z(o) && (o.$$moduleName = i), r.push([t, n, arguments]), c
                    }
                }
              }))
            }
          }))
        }(e))("ng", ["ngLocale"], ["$provide", function(e) {
          e.provider({
            $$sanitizeUri: Gr
          }), e.provider("$compile", vn).directive({
            a: Vi,
            input: wo,
            textarea: wo,
            form: zi,
            script: Ca,
            select: ka,
            option: Da,
            ngBind: Ao,
            ngBindHtml: ko,
            ngBindTemplate: To,
            ngClass: Eo,
            ngClassEven: No,
            ngClassOdd: Oo,
            ngCloak: Po,
            ngController: qo,
            ngForm: Yi,
            ngHide: ma,
            ngIf: jo,
            ngInclude: _o,
            ngInit: Uo,
            ngNonBindable: oa,
            ngPluralize: la,
            ngRef: da,
            ngRepeat: pa,
            ngShow: ga,
            ngStyle: va,
            ngSwitch: $a,
            ngSwitchWhen: ba,
            ngSwitchDefault: ya,
            ngOptions: ua,
            ngTransclude: xa,
            ngModel: ea,
            ngList: Ho,
            ngChange: Do,
            pattern: Ea,
            ngPattern: Ea,
            required: Ma,
            ngRequired: Ma,
            minlength: Na,
            ngMinlength: Na,
            maxlength: Oa,
            ngMaxlength: Oa,
            ngValue: So,
            ngModelOptions: ra
          }).directive({
            ngInclude: Vo,
            input: xo
          }).directive(Ui).directive(Io), e.provider({
            $anchorScroll: Zt,
            $animate: sn,
            $animateCss: cn,
            $$animateJs: on,
            $$animateQueue: an,
            $$AnimateRunner: ln,
            $$animateAsyncRun: un,
            $browser: pn,
            $cacheFactory: fn,
            $controller: kn,
            $document: Dn,
            $$isDocumentHidden: Mn,
            $exceptionHandler: En,
            $filter: mi,
            $$forceReflow: On,
            $interpolate: Qn,
            $interval: Jn,
            $$intervalFactory: Zn,
            $http: Wn,
            $httpParamSerializer: _n,
            $httpParamSerializerJQLike: Vn,
            $httpBackend: Yn,
            $xhrFactory: zn,
            $jsonpCallbacks: er,
            $location: vr,
            $log: $r,
            $parse: Lr,
            $rootScope: Br,
            $q: Rr,
            $$q: jr,
            $sce: Jr,
            $sceDelegate: Xr,
            $sniffer: Zr,
            $$taskTrackerFactory: ei,
            $templateCache: hn,
            $templateRequest: ri,
            $$testability: ii,
            $timeout: ai,
            $window: fi,
            $$rAF: Fr,
            $$jqLite: jt,
            $$Map: Ft,
            $$cookieReader: gi
          })
        }]).info({
          angularVersion: "1.8.2"
        })
      }(y), y.module("ngLocale", [], ["$provide", function(e) {
        var t = "one",
          n = "other";
        e.value("$locale", {
          DATETIME_FORMATS: {
            AMPMS: ["AM", "PM"],
            DAY: ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"],
            ERANAMES: ["Before Christ", "Anno Domini"],
            ERAS: ["BC", "AD"],
            FIRSTDAYOFWEEK: 6,
            MONTH: ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"],
            SHORTDAY: ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"],
            SHORTMONTH: ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"],
            STANDALONEMONTH: ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November",
              "December"
            ],
            WEEKENDRANGE: [5, 6],
            fullDate: "EEEE, MMMM d, y",
            longDate: "MMMM d, y",
            medium: "MMM d, y h:mm:ss a",
            mediumDate: "MMM d, y",
            mediumTime: "h:mm:ss a",
            short: "M/d/yy h:mm a",
            shortDate: "M/d/yy",
            shortTime: "h:mm a"
          },
          NUMBER_FORMATS: {
            CURRENCY_SYM: "$",
            DECIMAL_SEP: ".",
            GROUP_SEP: ",",
            PATTERNS: [{
              gSize: 3,
              lgSize: 3,
              maxFrac: 3,
              minFrac: 0,
              minInt: 1,
              negPre: "-",
              negSuf: "",
              posPre: "",
              posSuf: ""
            }, {
              gSize: 3,
              lgSize: 3,
              maxFrac: 2,
              minFrac: 2,
              minInt: 1,
              negPre: "-¤",
              negSuf: "",
              posPre: "¤",
              posSuf: ""
            }]
          },
          id: "en-us",
          localeID: "en_US",
          pluralCat: function(e, r) {
            var i = 0 | e,
              o = function(e, t) {
                var n = t;
                void 0 === n && (n = Math.min(function(e) {
                  var t = (e += "").indexOf(".");
                  return -1 == t ? 0 : e.length - t - 1
                }(e), 3));
                var r = Math.pow(10, n);
                return {
                  v: n,
                  f: (e * r | 0) % r
                }
              }(e, r);
            return 1 == i && 0 == o.v ? t : n
          }
        })
      }]), a((function() {
        De(e.document, Me)
      })))
    }(window), !window.angular.$$csp().noInlineStyle && window.angular.element(document.head).prepend(window.angular.element("<style>").text(
      '@charset "UTF-8";[ng\\:cloak],[ng-cloak],[data-ng-cloak],[x-ng-cloak],.ng-cloak,.x-ng-cloak,.ng-hide:not(.ng-hide-animate){display:none !important;}ng\\:form{display:block;}.ng-animate-shim{visibility:hidden;}.ng-anchor{position:absolute;}'
      ))
  }, {}],
  32: [function(e, t, n) {
    e("./angular"), t.exports = angular
  }, {
    "./angular": 31
  }],
  33: [function(e, t, n) {
    /*!
     * jQuery JavaScript Library v3.6.0
     * https://jquery.com/
     *
     * Includes Sizzle.js
     * https://sizzlejs.com/
     *
     * Copyright OpenJS Foundation and other contributors
     * Released under the MIT license
     * https://jquery.org/license
     *
     * Date: 2021-03-02T17:08Z
     */
    ! function(e, n) {
      "use strict";
      "object" == typeof t && "object" == typeof t.exports ? t.exports = e.document ? n(e, !0) : function(e) {
        if (!e.document) throw new Error("jQuery requires a window with a document");
        return n(e)
      } : n(e)
    }("undefined" != typeof window ? window : this, (function(e, t) {
      "use strict";
      var n = [],
        r = Object.getPrototypeOf,
        i = n.slice,
        o = n.flat ? function(e) {
          return n.flat.call(e)
        } : function(e) {
          return n.concat.apply([], e)
        },
        a = n.push,
        s = n.indexOf,
        u = {},
        l = u.toString,
        c = u.hasOwnProperty,
        d = c.toString,
        p = d.call(Object),
        f = {},
        h = function(e) {
          return "function" == typeof e && "number" != typeof e.nodeType && "function" != typeof e.item
        },
        g = function(e) {
          return null != e && e === e.window
        },
        m = e.document,
        v = {
          type: !0,
          src: !0,
          nonce: !0,
          noModule: !0
        };

      function $(e, t, n) {
        var r, i, o = (n = n || m).createElement("script");
        if (o.text = e, t)
          for (r in v)(i = t[r] || t.getAttribute && t.getAttribute(r)) && o.setAttribute(r, i);
        n.head.appendChild(o).parentNode.removeChild(o)
      }

      function b(e) {
        return null == e ? e + "" : "object" == typeof e || "function" == typeof e ? u[l.call(e)] || "object" : typeof e
      }
      var y = "3.6.0",
        w = function(e, t) {
          return new w.fn.init(e, t)
        };

      function x(e) {
        var t = !!e && "length" in e && e.length,
          n = b(e);
        return !h(e) && !g(e) && ("array" === n || 0 === t || "number" == typeof t && t > 0 && t - 1 in e)
      }
      w.fn = w.prototype = {
        jquery: y,
        constructor: w,
        length: 0,
        toArray: function() {
          return i.call(this)
        },
        get: function(e) {
          return null == e ? i.call(this) : e < 0 ? this[e + this.length] : this[e]
        },
        pushStack: function(e) {
          var t = w.merge(this.constructor(), e);
          return t.prevObject = this, t
        },
        each: function(e) {
          return w.each(this, e)
        },
        map: function(e) {
          return this.pushStack(w.map(this, (function(t, n) {
            return e.call(t, n, t)
          })))
        },
        slice: function() {
          return this.pushStack(i.apply(this, arguments))
        },
        first: function() {
          return this.eq(0)
        },
        last: function() {
          return this.eq(-1)
        },
        even: function() {
          return this.pushStack(w.grep(this, (function(e, t) {
            return (t + 1) % 2
          })))
        },
        odd: function() {
          return this.pushStack(w.grep(this, (function(e, t) {
            return t % 2
          })))
        },
        eq: function(e) {
          var t = this.length,
            n = +e + (e < 0 ? t : 0);
          return this.pushStack(n >= 0 && n < t ? [this[n]] : [])
        },
        end: function() {
          return this.prevObject || this.constructor()
        },
        push: a,
        sort: n.sort,
        splice: n.splice
      }, w.extend = w.fn.extend = function() {
        var e, t, n, r, i, o, a = arguments[0] || {},
          s = 1,
          u = arguments.length,
          l = !1;
        for ("boolean" == typeof a && (l = a, a = arguments[s] || {}, s++), "object" == typeof a || h(a) || (a = {}), s === u && (a = this,
          s--); s < u; s++)
          if (null != (e = arguments[s]))
            for (t in e) r = e[t], "__proto__" !== t && a !== r && (l && r && (w.isPlainObject(r) || (i = Array.isArray(r))) ? (n = a[t], o = i && !
              Array.isArray(n) ? [] : i || w.isPlainObject(n) ? n : {}, i = !1, a[t] = w.extend(l, o, r)) : void 0 !== r && (a[t] = r));
        return a
      }, w.extend({
        expando: "jQuery" + (y + Math.random()).replace(/\D/g, ""),
        isReady: !0,
        error: function(e) {
          throw new Error(e)
        },
        noop: function() {},
        isPlainObject: function(e) {
          var t, n;
          return !(!e || "[object Object]" !== l.call(e)) && (!(t = r(e)) || "function" == typeof(n = c.call(t, "constructor") && t
            .constructor) && d.call(n) === p)
        },
        isEmptyObject: function(e) {
          var t;
          for (t in e) return !1;
          return !0
        },
        globalEval: function(e, t, n) {
          $(e, {
            nonce: t && t.nonce
          }, n)
        },
        each: function(e, t) {
          var n, r = 0;
          if (x(e))
            for (n = e.length; r < n && !1 !== t.call(e[r], r, e[r]); r++);
          else
            for (r in e)
              if (!1 === t.call(e[r], r, e[r])) break;
          return e
        },
        makeArray: function(e, t) {
          var n = t || [];
          return null != e && (x(Object(e)) ? w.merge(n, "string" == typeof e ? [e] : e) : a.call(n, e)), n
        },
        inArray: function(e, t, n) {
          return null == t ? -1 : s.call(t, e, n)
        },
        merge: function(e, t) {
          for (var n = +t.length, r = 0, i = e.length; r < n; r++) e[i++] = t[r];
          return e.length = i, e
        },
        grep: function(e, t, n) {
          for (var r = [], i = 0, o = e.length, a = !n; i < o; i++) !t(e[i], i) !== a && r.push(e[i]);
          return r
        },
        map: function(e, t, n) {
          var r, i, a = 0,
            s = [];
          if (x(e))
            for (r = e.length; a < r; a++) null != (i = t(e[a], a, n)) && s.push(i);
          else
            for (a in e) null != (i = t(e[a], a, n)) && s.push(i);
          return o(s)
        },
        guid: 1,
        support: f
      }), "function" == typeof Symbol && (w.fn[Symbol.iterator] = n[Symbol.iterator]), w.each(
        "Boolean Number String Function Array Date RegExp Object Error Symbol".split(" "), (function(e, t) {
          u["[object " + t + "]"] = t.toLowerCase()
        }));
      var C =
        /*!
         * Sizzle CSS Selector Engine v2.3.6
         * https://sizzlejs.com/
         *
         * Copyright JS Foundation and other contributors
         * Released under the MIT license
         * https://js.foundation/
         *
         * Date: 2021-02-16
         */
        function(e) {
          var t, n, r, i, o, a, s, u, l, c, d, p, f, h, g, m, v, $, b, y = "sizzle" + 1 * new Date,
            w = e.document,
            x = 0,
            C = 0,
            S = ue(),
            A = ue(),
            T = ue(),
            k = ue(),
            D = function(e, t) {
              return e === t && (d = !0), 0
            },
            M = {}.hasOwnProperty,
            E = [],
            O = E.pop,
            N = E.push,
            P = E.push,
            q = E.slice,
            I = function(e, t) {
              for (var n = 0, r = e.length; n < r; n++)
                if (e[n] === t) return n;
              return -1
            },
            L = "checked|selected|async|autofocus|autoplay|controls|defer|disabled|hidden|ismap|loop|multiple|open|readonly|required|scoped",
            R = "[\\x20\\t\\r\\n\\f]",
            j = "(?:\\\\[\\da-fA-F]{1,6}[\\x20\\t\\r\\n\\f]?|\\\\[^\\r\\n\\f]|[\\w-]|[^\0-\\x7f])+",
            _ = "\\[[\\x20\\t\\r\\n\\f]*(" + j + ")(?:" + R + "*([*^$|!~]?=)" + R + "*(?:'((?:\\\\.|[^\\\\'])*)'|\"((?:\\\\.|[^\\\\\"])*)\"|(" + j +
            "))|)" + R + "*\\]",
            V = ":(" + j + ")(?:\\((('((?:\\\\.|[^\\\\'])*)'|\"((?:\\\\.|[^\\\\\"])*)\")|((?:\\\\.|[^\\\\()[\\]]|" + _ + ")*)|.*)\\)|)",
            U = new RegExp(R + "+", "g"),
            H = new RegExp("^[\\x20\\t\\r\\n\\f]+|((?:^|[^\\\\])(?:\\\\.)*)[\\x20\\t\\r\\n\\f]+$", "g"),
            F = new RegExp("^[\\x20\\t\\r\\n\\f]*,[\\x20\\t\\r\\n\\f]*"),
            B = new RegExp("^[\\x20\\t\\r\\n\\f]*([>+~]|[\\x20\\t\\r\\n\\f])[\\x20\\t\\r\\n\\f]*"),
            G = new RegExp(R + "|>"),
            W = new RegExp(V),
            z = new RegExp("^" + j + "$"),
            Y = {
              ID: new RegExp("^#(" + j + ")"),
              CLASS: new RegExp("^\\.(" + j + ")"),
              TAG: new RegExp("^(" + j + "|[*])"),
              ATTR: new RegExp("^" + _),
              PSEUDO: new RegExp("^" + V),
              CHILD: new RegExp(
                "^:(only|first|last|nth|nth-last)-(child|of-type)(?:\\([\\x20\\t\\r\\n\\f]*(even|odd|(([+-]|)(\\d*)n|)[\\x20\\t\\r\\n\\f]*(?:([+-]|)[\\x20\\t\\r\\n\\f]*(\\d+)|))[\\x20\\t\\r\\n\\f]*\\)|)",
                "i"),
              bool: new RegExp("^(?:" + L + ")$", "i"),
              needsContext: new RegExp(
                "^[\\x20\\t\\r\\n\\f]*[>+~]|:(even|odd|eq|gt|lt|nth|first|last)(?:\\([\\x20\\t\\r\\n\\f]*((?:-\\d)?\\d*)[\\x20\\t\\r\\n\\f]*\\)|)(?=[^-]|$)",
                "i")
            },
            K = /HTML$/i,
            Q = /^(?:input|select|textarea|button)$/i,
            X = /^h\d$/i,
            J = /^[^{]+\{\s*\[native \w/,
            Z = /^(?:#([\w-]+)|(\w+)|\.([\w-]+))$/,
            ee = /[+~]/,
            te = new RegExp("\\\\[\\da-fA-F]{1,6}[\\x20\\t\\r\\n\\f]?|\\\\([^\\r\\n\\f])", "g"),
            ne = function(e, t) {
              var n = "0x" + e.slice(1) - 65536;
              return t || (n < 0 ? String.fromCharCode(n + 65536) : String.fromCharCode(n >> 10 | 55296, 1023 & n | 56320))
            },
            re = /([\0-\x1f\x7f]|^-?\d)|^-$|[^\0-\x1f\x7f-\uFFFF\w-]/g,
            ie = function(e, t) {
              return t ? "\0" === e ? "�" : e.slice(0, -1) + "\\" + e.charCodeAt(e.length - 1).toString(16) + " " : "\\" + e
            },
            oe = function() {
              p()
            },
            ae = ye((function(e) {
              return !0 === e.disabled && "fieldset" === e.nodeName.toLowerCase()
            }), {
              dir: "parentNode",
              next: "legend"
            });
          try {
            P.apply(E = q.call(w.childNodes), w.childNodes), E[w.childNodes.length].nodeType
          } catch (e) {
            P = {
              apply: E.length ? function(e, t) {
                N.apply(e, q.call(t))
              } : function(e, t) {
                for (var n = e.length, r = 0; e[n++] = t[r++];);
                e.length = n - 1
              }
            }
          }

          function se(e, t, r, i) {
            var o, s, l, c, d, h, v, $ = t && t.ownerDocument,
              w = t ? t.nodeType : 9;
            if (r = r || [], "string" != typeof e || !e || 1 !== w && 9 !== w && 11 !== w) return r;
            if (!i && (p(t), t = t || f, g)) {
              if (11 !== w && (d = Z.exec(e)))
                if (o = d[1]) {
                  if (9 === w) {
                    if (!(l = t.getElementById(o))) return r;
                    if (l.id === o) return r.push(l), r
                  } else if ($ && (l = $.getElementById(o)) && b(t, l) && l.id === o) return r.push(l), r
                } else {
                  if (d[2]) return P.apply(r, t.getElementsByTagName(e)), r;
                  if ((o = d[3]) && n.getElementsByClassName && t.getElementsByClassName) return P.apply(r, t.getElementsByClassName(o)), r
                } if (n.qsa && !k[e + " "] && (!m || !m.test(e)) && (1 !== w || "object" !== t.nodeName.toLowerCase())) {
                if (v = e, $ = t, 1 === w && (G.test(e) || B.test(e))) {
                  for (($ = ee.test(e) && ve(t.parentNode) || t) === t && n.scope || ((c = t.getAttribute("id")) ? c = c.replace(re, ie) : t
                      .setAttribute("id", c = y)), s = (h = a(e)).length; s--;) h[s] = (c ? "#" + c : ":scope") + " " + be(h[s]);
                  v = h.join(",")
                }
                try {
                  return P.apply(r, $.querySelectorAll(v)), r
                } catch (t) {
                  k(e, !0)
                } finally {
                  c === y && t.removeAttribute("id")
                }
              }
            }
            return u(e.replace(H, "$1"), t, r, i)
          }

          function ue() {
            var e = [];
            return function t(n, i) {
              return e.push(n + " ") > r.cacheLength && delete t[e.shift()], t[n + " "] = i
            }
          }

          function le(e) {
            return e[y] = !0, e
          }

          function ce(e) {
            var t = f.createElement("fieldset");
            try {
              return !!e(t)
            } catch (e) {
              return !1
            } finally {
              t.parentNode && t.parentNode.removeChild(t), t = null
            }
          }

          function de(e, t) {
            for (var n = e.split("|"), i = n.length; i--;) r.attrHandle[n[i]] = t
          }

          function pe(e, t) {
            var n = t && e,
              r = n && 1 === e.nodeType && 1 === t.nodeType && e.sourceIndex - t.sourceIndex;
            if (r) return r;
            if (n)
              for (; n = n.nextSibling;)
                if (n === t) return -1;
            return e ? 1 : -1
          }

          function fe(e) {
            return function(t) {
              return "input" === t.nodeName.toLowerCase() && t.type === e
            }
          }

          function he(e) {
            return function(t) {
              var n = t.nodeName.toLowerCase();
              return ("input" === n || "button" === n) && t.type === e
            }
          }

          function ge(e) {
            return function(t) {
              return "form" in t ? t.parentNode && !1 === t.disabled ? "label" in t ? "label" in t.parentNode ? t.parentNode.disabled === e : t
                .disabled === e : t.isDisabled === e || t.isDisabled !== !e && ae(t) === e : t.disabled === e : "label" in t && t.disabled === e
            }
          }

          function me(e) {
            return le((function(t) {
              return t = +t, le((function(n, r) {
                for (var i, o = e([], n.length, t), a = o.length; a--;) n[i = o[a]] && (n[i] = !(r[i] = n[i]))
              }))
            }))
          }

          function ve(e) {
            return e && void 0 !== e.getElementsByTagName && e
          }
          for (t in n = se.support = {}, o = se.isXML = function(e) {
              var t = e && e.namespaceURI,
                n = e && (e.ownerDocument || e).documentElement;
              return !K.test(t || n && n.nodeName || "HTML")
            }, p = se.setDocument = function(e) {
              var t, i, a = e ? e.ownerDocument || e : w;
              return a != f && 9 === a.nodeType && a.documentElement ? (h = (f = a).documentElement, g = !o(f), w != f && (i = f.defaultView) && i
                .top !== i && (i.addEventListener ? i.addEventListener("unload", oe, !1) : i.attachEvent && i.attachEvent("onunload", oe)), n
                .scope = ce((function(e) {
                  return h.appendChild(e).appendChild(f.createElement("div")), void 0 !== e.querySelectorAll && !e.querySelectorAll(
                    ":scope fieldset div").length
                })), n.attributes = ce((function(e) {
                  return e.className = "i", !e.getAttribute("className")
                })), n.getElementsByTagName = ce((function(e) {
                  return e.appendChild(f.createComment("")), !e.getElementsByTagName("*").length
                })), n.getElementsByClassName = J.test(f.getElementsByClassName), n.getById = ce((function(e) {
                  return h.appendChild(e).id = y, !f.getElementsByName || !f.getElementsByName(y).length
                })), n.getById ? (r.filter.ID = function(e) {
                  var t = e.replace(te, ne);
                  return function(e) {
                    return e.getAttribute("id") === t
                  }
                }, r.find.ID = function(e, t) {
                  if (void 0 !== t.getElementById && g) {
                    var n = t.getElementById(e);
                    return n ? [n] : []
                  }
                }) : (r.filter.ID = function(e) {
                  var t = e.replace(te, ne);
                  return function(e) {
                    var n = void 0 !== e.getAttributeNode && e.getAttributeNode("id");
                    return n && n.value === t
                  }
                }, r.find.ID = function(e, t) {
                  if (void 0 !== t.getElementById && g) {
                    var n, r, i, o = t.getElementById(e);
                    if (o) {
                      if ((n = o.getAttributeNode("id")) && n.value === e) return [o];
                      for (i = t.getElementsByName(e), r = 0; o = i[r++];)
                        if ((n = o.getAttributeNode("id")) && n.value === e) return [o]
                    }
                    return []
                  }
                }), r.find.TAG = n.getElementsByTagName ? function(e, t) {
                  return void 0 !== t.getElementsByTagName ? t.getElementsByTagName(e) : n.qsa ? t.querySelectorAll(e) : void 0
                } : function(e, t) {
                  var n, r = [],
                    i = 0,
                    o = t.getElementsByTagName(e);
                  if ("*" === e) {
                    for (; n = o[i++];) 1 === n.nodeType && r.push(n);
                    return r
                  }
                  return o
                }, r.find.CLASS = n.getElementsByClassName && function(e, t) {
                  if (void 0 !== t.getElementsByClassName && g) return t.getElementsByClassName(e)
                }, v = [], m = [], (n.qsa = J.test(f.querySelectorAll)) && (ce((function(e) {
                  var t;
                  h.appendChild(e).innerHTML = "<a id='" + y + "'></a><select id='" + y +
                    "-\r\\' msallowcapture=''><option selected=''></option></select>", e.querySelectorAll("[msallowcapture^='']").length && m
                    .push("[*^$]=[\\x20\\t\\r\\n\\f]*(?:''|\"\")"), e.querySelectorAll("[selected]").length || m.push(
                      "\\[[\\x20\\t\\r\\n\\f]*(?:value|" + L + ")"), e.querySelectorAll("[id~=" + y + "-]").length || m.push("~="), (t = f
                      .createElement("input")).setAttribute("name", ""), e.appendChild(t), e.querySelectorAll("[name='']").length || m.push(
                      "\\[[\\x20\\t\\r\\n\\f]*name[\\x20\\t\\r\\n\\f]*=[\\x20\\t\\r\\n\\f]*(?:''|\"\")"), e.querySelectorAll(":checked")
                    .length || m.push(":checked"), e.querySelectorAll("a#" + y + "+*").length || m.push(".#.+[+~]"), e.querySelectorAll(
                      "\\\f"), m.push("[\\r\\n\\f]")
                })), ce((function(e) {
                  e.innerHTML = "<a href='' disabled='disabled'></a><select disabled='disabled'><option/></select>";
                  var t = f.createElement("input");
                  t.setAttribute("type", "hidden"), e.appendChild(t).setAttribute("name", "D"), e.querySelectorAll("[name=d]").length && m
                    .push("name[\\x20\\t\\r\\n\\f]*[*^$|!~]?="), 2 !== e.querySelectorAll(":enabled").length && m.push(":enabled",
                      ":disabled"), h.appendChild(e).disabled = !0, 2 !== e.querySelectorAll(":disabled").length && m.push(":enabled",
                      ":disabled"), e.querySelectorAll("*,:x"), m.push(",.*:")
                }))), (n.matchesSelector = J.test($ = h.matches || h.webkitMatchesSelector || h.mozMatchesSelector || h.oMatchesSelector || h
                  .msMatchesSelector)) && ce((function(e) {
                  n.disconnectedMatch = $.call(e, "*"), $.call(e, "[s!='']:x"), v.push("!=", V)
                })), m = m.length && new RegExp(m.join("|")), v = v.length && new RegExp(v.join("|")), t = J.test(h.compareDocumentPosition), b =
                t || J.test(h.contains) ? function(e, t) {
                  var n = 9 === e.nodeType ? e.documentElement : e,
                    r = t && t.parentNode;
                  return e === r || !(!r || 1 !== r.nodeType || !(n.contains ? n.contains(r) : e.compareDocumentPosition && 16 & e
                    .compareDocumentPosition(r)))
                } : function(e, t) {
                  if (t)
                    for (; t = t.parentNode;)
                      if (t === e) return !0;
                  return !1
                }, D = t ? function(e, t) {
                  if (e === t) return d = !0, 0;
                  var r = !e.compareDocumentPosition - !t.compareDocumentPosition;
                  return r || (1 & (r = (e.ownerDocument || e) == (t.ownerDocument || t) ? e.compareDocumentPosition(t) : 1) || !n.sortDetached && t
                    .compareDocumentPosition(e) === r ? e == f || e.ownerDocument == w && b(w, e) ? -1 : t == f || t.ownerDocument == w && b(w,
                    t) ? 1 : c ? I(c, e) - I(c, t) : 0 : 4 & r ? -1 : 1)
                } : function(e, t) {
                  if (e === t) return d = !0, 0;
                  var n, r = 0,
                    i = e.parentNode,
                    o = t.parentNode,
                    a = [e],
                    s = [t];
                  if (!i || !o) return e == f ? -1 : t == f ? 1 : i ? -1 : o ? 1 : c ? I(c, e) - I(c, t) : 0;
                  if (i === o) return pe(e, t);
                  for (n = e; n = n.parentNode;) a.unshift(n);
                  for (n = t; n = n.parentNode;) s.unshift(n);
                  for (; a[r] === s[r];) r++;
                  return r ? pe(a[r], s[r]) : a[r] == w ? -1 : s[r] == w ? 1 : 0
                }, f) : f
            }, se.matches = function(e, t) {
              return se(e, null, null, t)
            }, se.matchesSelector = function(e, t) {
              if (p(e), n.matchesSelector && g && !k[t + " "] && (!v || !v.test(t)) && (!m || !m.test(t))) try {
                var r = $.call(e, t);
                if (r || n.disconnectedMatch || e.document && 11 !== e.document.nodeType) return r
              } catch (e) {
                k(t, !0)
              }
              return se(t, f, null, [e]).length > 0
            }, se.contains = function(e, t) {
              return (e.ownerDocument || e) != f && p(e), b(e, t)
            }, se.attr = function(e, t) {
              (e.ownerDocument || e) != f && p(e);
              var i = r.attrHandle[t.toLowerCase()],
                o = i && M.call(r.attrHandle, t.toLowerCase()) ? i(e, t, !g) : void 0;
              return void 0 !== o ? o : n.attributes || !g ? e.getAttribute(t) : (o = e.getAttributeNode(t)) && o.specified ? o.value : null
            }, se.escape = function(e) {
              return (e + "").replace(re, ie)
            }, se.error = function(e) {
              throw new Error("Syntax error, unrecognized expression: " + e)
            }, se.uniqueSort = function(e) {
              var t, r = [],
                i = 0,
                o = 0;
              if (d = !n.detectDuplicates, c = !n.sortStable && e.slice(0), e.sort(D), d) {
                for (; t = e[o++];) t === e[o] && (i = r.push(o));
                for (; i--;) e.splice(r[i], 1)
              }
              return c = null, e
            }, i = se.getText = function(e) {
              var t, n = "",
                r = 0,
                o = e.nodeType;
              if (o) {
                if (1 === o || 9 === o || 11 === o) {
                  if ("string" == typeof e.textContent) return e.textContent;
                  for (e = e.firstChild; e; e = e.nextSibling) n += i(e)
                } else if (3 === o || 4 === o) return e.nodeValue
              } else
                for (; t = e[r++];) n += i(t);
              return n
            }, (r = se.selectors = {
              cacheLength: 50,
              createPseudo: le,
              match: Y,
              attrHandle: {},
              find: {},
              relative: {
                ">": {
                  dir: "parentNode",
                  first: !0
                },
                " ": {
                  dir: "parentNode"
                },
                "+": {
                  dir: "previousSibling",
                  first: !0
                },
                "~": {
                  dir: "previousSibling"
                }
              },
              preFilter: {
                ATTR: function(e) {
                  return e[1] = e[1].replace(te, ne), e[3] = (e[3] || e[4] || e[5] || "").replace(te, ne), "~=" === e[2] && (e[3] = " " + e[3] +
                    " "), e.slice(0, 4)
                },
                CHILD: function(e) {
                  return e[1] = e[1].toLowerCase(), "nth" === e[1].slice(0, 3) ? (e[3] || se.error(e[0]), e[4] = +(e[4] ? e[5] + (e[6] || 1) : 2 *
                    ("even" === e[3] || "odd" === e[3])), e[5] = +(e[7] + e[8] || "odd" === e[3])) : e[3] && se.error(e[0]), e
                },
                PSEUDO: function(e) {
                  var t, n = !e[6] && e[2];
                  return Y.CHILD.test(e[0]) ? null : (e[3] ? e[2] = e[4] || e[5] || "" : n && W.test(n) && (t = a(n, !0)) && (t = n.indexOf(")", n
                    .length - t) - n.length) && (e[0] = e[0].slice(0, t), e[2] = n.slice(0, t)), e.slice(0, 3))
                }
              },
              filter: {
                TAG: function(e) {
                  var t = e.replace(te, ne).toLowerCase();
                  return "*" === e ? function() {
                    return !0
                  } : function(e) {
                    return e.nodeName && e.nodeName.toLowerCase() === t
                  }
                },
                CLASS: function(e) {
                  var t = S[e + " "];
                  return t || (t = new RegExp("(^|[\\x20\\t\\r\\n\\f])" + e + "(" + R + "|$)")) && S(e, (function(e) {
                    return t.test("string" == typeof e.className && e.className || void 0 !== e.getAttribute && e.getAttribute("class") ||
                      "")
                  }))
                },
                ATTR: function(e, t, n) {
                  return function(r) {
                    var i = se.attr(r, e);
                    return null == i ? "!=" === t : !t || (i += "", "=" === t ? i === n : "!=" === t ? i !== n : "^=" === t ? n && 0 === i
                      .indexOf(n) : "*=" === t ? n && i.indexOf(n) > -1 : "$=" === t ? n && i.slice(-n.length) === n : "~=" === t ? (" " + i
                        .replace(U, " ") + " ").indexOf(n) > -1 : "|=" === t && (i === n || i.slice(0, n.length + 1) === n + "-"))
                  }
                },
                CHILD: function(e, t, n, r, i) {
                  var o = "nth" !== e.slice(0, 3),
                    a = "last" !== e.slice(-4),
                    s = "of-type" === t;
                  return 1 === r && 0 === i ? function(e) {
                    return !!e.parentNode
                  } : function(t, n, u) {
                    var l, c, d, p, f, h, g = o !== a ? "nextSibling" : "previousSibling",
                      m = t.parentNode,
                      v = s && t.nodeName.toLowerCase(),
                      $ = !u && !s,
                      b = !1;
                    if (m) {
                      if (o) {
                        for (; g;) {
                          for (p = t; p = p[g];)
                            if (s ? p.nodeName.toLowerCase() === v : 1 === p.nodeType) return !1;
                          h = g = "only" === e && !h && "nextSibling"
                        }
                        return !0
                      }
                      if (h = [a ? m.firstChild : m.lastChild], a && $) {
                        for (b = (f = (l = (c = (d = (p = m)[y] || (p[y] = {}))[p.uniqueID] || (d[p.uniqueID] = {}))[e] || [])[0] === x && l[
                          1]) && l[2], p = f && m.childNodes[f]; p = ++f && p && p[g] || (b = f = 0) || h.pop();)
                          if (1 === p.nodeType && ++b && p === t) {
                            c[e] = [x, f, b];
                            break
                          }
                      } else if ($ && (b = f = (l = (c = (d = (p = t)[y] || (p[y] = {}))[p.uniqueID] || (d[p.uniqueID] = {}))[e] || [])[0] ===
                          x && l[1]), !1 === b)
                        for (;
                          (p = ++f && p && p[g] || (b = f = 0) || h.pop()) && ((s ? p.nodeName.toLowerCase() !== v : 1 !== p.nodeType) || !++
                            b || ($ && ((c = (d = p[y] || (p[y] = {}))[p.uniqueID] || (d[p.uniqueID] = {}))[e] = [x, b]), p !== t)););
                      return (b -= i) === r || b % r == 0 && b / r >= 0
                    }
                  }
                },
                PSEUDO: function(e, t) {
                  var n, i = r.pseudos[e] || r.setFilters[e.toLowerCase()] || se.error("unsupported pseudo: " + e);
                  return i[y] ? i(t) : i.length > 1 ? (n = [e, e, "", t], r.setFilters.hasOwnProperty(e.toLowerCase()) ? le((function(e, n) {
                    for (var r, o = i(e, t), a = o.length; a--;) e[r = I(e, o[a])] = !(n[r] = o[a])
                  })) : function(e) {
                    return i(e, 0, n)
                  }) : i
                }
              },
              pseudos: {
                not: le((function(e) {
                  var t = [],
                    n = [],
                    r = s(e.replace(H, "$1"));
                  return r[y] ? le((function(e, t, n, i) {
                    for (var o, a = r(e, null, i, []), s = e.length; s--;)(o = a[s]) && (e[s] = !(t[s] = o))
                  })) : function(e, i, o) {
                    return t[0] = e, r(t, null, o, n), t[0] = null, !n.pop()
                  }
                })),
                has: le((function(e) {
                  return function(t) {
                    return se(e, t).length > 0
                  }
                })),
                contains: le((function(e) {
                  return e = e.replace(te, ne),
                    function(t) {
                      return (t.textContent || i(t)).indexOf(e) > -1
                    }
                })),
                lang: le((function(e) {
                  return z.test(e || "") || se.error("unsupported lang: " + e), e = e.replace(te, ne).toLowerCase(),
                    function(t) {
                      var n;
                      do {
                        if (n = g ? t.lang : t.getAttribute("xml:lang") || t.getAttribute("lang")) return (n = n.toLowerCase()) === e || 0 ===
                          n.indexOf(e + "-")
                      } while ((t = t.parentNode) && 1 === t.nodeType);
                      return !1
                    }
                })),
                target: function(t) {
                  var n = e.location && e.location.hash;
                  return n && n.slice(1) === t.id
                },
                root: function(e) {
                  return e === h
                },
                focus: function(e) {
                  return e === f.activeElement && (!f.hasFocus || f.hasFocus()) && !!(e.type || e.href || ~e.tabIndex)
                },
                enabled: ge(!1),
                disabled: ge(!0),
                checked: function(e) {
                  var t = e.nodeName.toLowerCase();
                  return "input" === t && !!e.checked || "option" === t && !!e.selected
                },
                selected: function(e) {
                  return e.parentNode && e.parentNode.selectedIndex, !0 === e.selected
                },
                empty: function(e) {
                  for (e = e.firstChild; e; e = e.nextSibling)
                    if (e.nodeType < 6) return !1;
                  return !0
                },
                parent: function(e) {
                  return !r.pseudos.empty(e)
                },
                header: function(e) {
                  return X.test(e.nodeName)
                },
                input: function(e) {
                  return Q.test(e.nodeName)
                },
                button: function(e) {
                  var t = e.nodeName.toLowerCase();
                  return "input" === t && "button" === e.type || "button" === t
                },
                text: function(e) {
                  var t;
                  return "input" === e.nodeName.toLowerCase() && "text" === e.type && (null == (t = e.getAttribute("type")) || "text" === t
                    .toLowerCase())
                },
                first: me((function() {
                  return [0]
                })),
                last: me((function(e, t) {
                  return [t - 1]
                })),
                eq: me((function(e, t, n) {
                  return [n < 0 ? n + t : n]
                })),
                even: me((function(e, t) {
                  for (var n = 0; n < t; n += 2) e.push(n);
                  return e
                })),
                odd: me((function(e, t) {
                  for (var n = 1; n < t; n += 2) e.push(n);
                  return e
                })),
                lt: me((function(e, t, n) {
                  for (var r = n < 0 ? n + t : n > t ? t : n; --r >= 0;) e.push(r);
                  return e
                })),
                gt: me((function(e, t, n) {
                  for (var r = n < 0 ? n + t : n; ++r < t;) e.push(r);
                  return e
                }))
              }
            }).pseudos.nth = r.pseudos.eq, {
              radio: !0,
              checkbox: !0,
              file: !0,
              password: !0,
              image: !0
            }) r.pseudos[t] = fe(t);
          for (t in {
              submit: !0,
              reset: !0
            }) r.pseudos[t] = he(t);

          function $e() {}

          function be(e) {
            for (var t = 0, n = e.length, r = ""; t < n; t++) r += e[t].value;
            return r
          }

          function ye(e, t, n) {
            var r = t.dir,
              i = t.next,
              o = i || r,
              a = n && "parentNode" === o,
              s = C++;
            return t.first ? function(t, n, i) {
              for (; t = t[r];)
                if (1 === t.nodeType || a) return e(t, n, i);
              return !1
            } : function(t, n, u) {
              var l, c, d, p = [x, s];
              if (u) {
                for (; t = t[r];)
                  if ((1 === t.nodeType || a) && e(t, n, u)) return !0
              } else
                for (; t = t[r];)
                  if (1 === t.nodeType || a)
                    if (c = (d = t[y] || (t[y] = {}))[t.uniqueID] || (d[t.uniqueID] = {}), i && i === t.nodeName.toLowerCase()) t = t[r] || t;
                    else {
                      if ((l = c[o]) && l[0] === x && l[1] === s) return p[2] = l[2];
                      if (c[o] = p, p[2] = e(t, n, u)) return !0
                    } return !1
            }
          }

          function we(e) {
            return e.length > 1 ? function(t, n, r) {
              for (var i = e.length; i--;)
                if (!e[i](t, n, r)) return !1;
              return !0
            } : e[0]
          }

          function xe(e, t, n, r, i) {
            for (var o, a = [], s = 0, u = e.length, l = null != t; s < u; s++)(o = e[s]) && (n && !n(o, r, i) || (a.push(o), l && t.push(s)));
            return a
          }

          function Ce(e, t, n, r, i, o) {
            return r && !r[y] && (r = Ce(r)), i && !i[y] && (i = Ce(i, o)), le((function(o, a, s, u) {
              var l, c, d, p = [],
                f = [],
                h = a.length,
                g = o || function(e, t, n) {
                  for (var r = 0, i = t.length; r < i; r++) se(e, t[r], n);
                  return n
                }(t || "*", s.nodeType ? [s] : s, []),
                m = !e || !o && t ? g : xe(g, p, e, s, u),
                v = n ? i || (o ? e : h || r) ? [] : a : m;
              if (n && n(m, v, s, u), r)
                for (l = xe(v, f), r(l, [], s, u), c = l.length; c--;)(d = l[c]) && (v[f[c]] = !(m[f[c]] = d));
              if (o) {
                if (i || e) {
                  if (i) {
                    for (l = [], c = v.length; c--;)(d = v[c]) && l.push(m[c] = d);
                    i(null, v = [], l, u)
                  }
                  for (c = v.length; c--;)(d = v[c]) && (l = i ? I(o, d) : p[c]) > -1 && (o[l] = !(a[l] = d))
                }
              } else v = xe(v === a ? v.splice(h, v.length) : v), i ? i(null, a, v, u) : P.apply(a, v)
            }))
          }

          function Se(e) {
            for (var t, n, i, o = e.length, a = r.relative[e[0].type], s = a || r.relative[" "], u = a ? 1 : 0, c = ye((function(e) {
                return e === t
              }), s, !0), d = ye((function(e) {
                return I(t, e) > -1
              }), s, !0), p = [function(e, n, r) {
                var i = !a && (r || n !== l) || ((t = n).nodeType ? c(e, n, r) : d(e, n, r));
                return t = null, i
              }]; u < o; u++)
              if (n = r.relative[e[u].type]) p = [ye(we(p), n)];
              else {
                if ((n = r.filter[e[u].type].apply(null, e[u].matches))[y]) {
                  for (i = ++u; i < o && !r.relative[e[i].type]; i++);
                  return Ce(u > 1 && we(p), u > 1 && be(e.slice(0, u - 1).concat({
                    value: " " === e[u - 2].type ? "*" : ""
                  })).replace(H, "$1"), n, u < i && Se(e.slice(u, i)), i < o && Se(e = e.slice(i)), i < o && be(e))
                }
                p.push(n)
              } return we(p)
          }
          return $e.prototype = r.filters = r.pseudos, r.setFilters = new $e, a = se.tokenize = function(e, t) {
            var n, i, o, a, s, u, l, c = A[e + " "];
            if (c) return t ? 0 : c.slice(0);
            for (s = e, u = [], l = r.preFilter; s;) {
              for (a in n && !(i = F.exec(s)) || (i && (s = s.slice(i[0].length) || s), u.push(o = [])), n = !1, (i = B.exec(s)) && (n = i.shift(),
                  o.push({
                    value: n,
                    type: i[0].replace(H, " ")
                  }), s = s.slice(n.length)), r.filter) !(i = Y[a].exec(s)) || l[a] && !(i = l[a](i)) || (n = i.shift(), o.push({
                value: n,
                type: a,
                matches: i
              }), s = s.slice(n.length));
              if (!n) break
            }
            return t ? s.length : s ? se.error(e) : A(e, u).slice(0)
          }, s = se.compile = function(e, t) {
            var n, i = [],
              o = [],
              s = T[e + " "];
            if (!s) {
              for (t || (t = a(e)), n = t.length; n--;)(s = Se(t[n]))[y] ? i.push(s) : o.push(s);
              (s = T(e, function(e, t) {
                var n = t.length > 0,
                  i = e.length > 0,
                  o = function(o, a, s, u, c) {
                    var d, h, m, v = 0,
                      $ = "0",
                      b = o && [],
                      y = [],
                      w = l,
                      C = o || i && r.find.TAG("*", c),
                      S = x += null == w ? 1 : Math.random() || .1,
                      A = C.length;
                    for (c && (l = a == f || a || c); $ !== A && null != (d = C[$]); $++) {
                      if (i && d) {
                        for (h = 0, a || d.ownerDocument == f || (p(d), s = !g); m = e[h++];)
                          if (m(d, a || f, s)) {
                            u.push(d);
                            break
                          } c && (x = S)
                      }
                      n && ((d = !m && d) && v--, o && b.push(d))
                    }
                    if (v += $, n && $ !== v) {
                      for (h = 0; m = t[h++];) m(b, y, a, s);
                      if (o) {
                        if (v > 0)
                          for (; $--;) b[$] || y[$] || (y[$] = O.call(u));
                        y = xe(y)
                      }
                      P.apply(u, y), c && !o && y.length > 0 && v + t.length > 1 && se.uniqueSort(u)
                    }
                    return c && (x = S, l = w), b
                  };
                return n ? le(o) : o
              }(o, i))).selector = e
            }
            return s
          }, u = se.select = function(e, t, n, i) {
            var o, u, l, c, d, p = "function" == typeof e && e,
              f = !i && a(e = p.selector || e);
            if (n = n || [], 1 === f.length) {
              if ((u = f[0] = f[0].slice(0)).length > 2 && "ID" === (l = u[0]).type && 9 === t.nodeType && g && r.relative[u[1].type]) {
                if (!(t = (r.find.ID(l.matches[0].replace(te, ne), t) || [])[0])) return n;
                p && (t = t.parentNode), e = e.slice(u.shift().value.length)
              }
              for (o = Y.needsContext.test(e) ? 0 : u.length; o-- && (l = u[o], !r.relative[c = l.type]);)
                if ((d = r.find[c]) && (i = d(l.matches[0].replace(te, ne), ee.test(u[0].type) && ve(t.parentNode) || t))) {
                  if (u.splice(o, 1), !(e = i.length && be(u))) return P.apply(n, i), n;
                  break
                }
            }
            return (p || s(e, f))(i, t, !g, n, !t || ee.test(e) && ve(t.parentNode) || t), n
          }, n.sortStable = y.split("").sort(D).join("") === y, n.detectDuplicates = !!d, p(), n.sortDetached = ce((function(e) {
            return 1 & e.compareDocumentPosition(f.createElement("fieldset"))
          })), ce((function(e) {
            return e.innerHTML = "<a href='#'></a>", "#" === e.firstChild.getAttribute("href")
          })) || de("type|href|height|width", (function(e, t, n) {
            if (!n) return e.getAttribute(t, "type" === t.toLowerCase() ? 1 : 2)
          })), n.attributes && ce((function(e) {
            return e.innerHTML = "<input/>", e.firstChild.setAttribute("value", ""), "" === e.firstChild.getAttribute("value")
          })) || de("value", (function(e, t, n) {
            if (!n && "input" === e.nodeName.toLowerCase()) return e.defaultValue
          })), ce((function(e) {
            return null == e.getAttribute("disabled")
          })) || de(L, (function(e, t, n) {
            var r;
            if (!n) return !0 === e[t] ? t.toLowerCase() : (r = e.getAttributeNode(t)) && r.specified ? r.value : null
          })), se
        }(e);
      w.find = C, w.expr = C.selectors, w.expr[":"] = w.expr.pseudos, w.uniqueSort = w.unique = C.uniqueSort, w.text = C.getText, w.isXMLDoc = C
        .isXML, w.contains = C.contains, w.escapeSelector = C.escape;
      var S = function(e, t, n) {
          for (var r = [], i = void 0 !== n;
            (e = e[t]) && 9 !== e.nodeType;)
            if (1 === e.nodeType) {
              if (i && w(e).is(n)) break;
              r.push(e)
            } return r
        },
        A = function(e, t) {
          for (var n = []; e; e = e.nextSibling) 1 === e.nodeType && e !== t && n.push(e);
          return n
        },
        T = w.expr.match.needsContext;

      function k(e, t) {
        return e.nodeName && e.nodeName.toLowerCase() === t.toLowerCase()
      }
      var D = /^<([a-z][^\/\0>:\x20\t\r\n\f]*)[\x20\t\r\n\f]*\/?>(?:<\/\1>|)$/i;

      function M(e, t, n) {
        return h(t) ? w.grep(e, (function(e, r) {
          return !!t.call(e, r, e) !== n
        })) : t.nodeType ? w.grep(e, (function(e) {
          return e === t !== n
        })) : "string" != typeof t ? w.grep(e, (function(e) {
          return s.call(t, e) > -1 !== n
        })) : w.filter(t, e, n)
      }
      w.filter = function(e, t, n) {
        var r = t[0];
        return n && (e = ":not(" + e + ")"), 1 === t.length && 1 === r.nodeType ? w.find.matchesSelector(r, e) ? [r] : [] : w.find.matches(e, w
          .grep(t, (function(e) {
            return 1 === e.nodeType
          })))
      }, w.fn.extend({
        find: function(e) {
          var t, n, r = this.length,
            i = this;
          if ("string" != typeof e) return this.pushStack(w(e).filter((function() {
            for (t = 0; t < r; t++)
              if (w.contains(i[t], this)) return !0
          })));
          for (n = this.pushStack([]), t = 0; t < r; t++) w.find(e, i[t], n);
          return r > 1 ? w.uniqueSort(n) : n
        },
        filter: function(e) {
          return this.pushStack(M(this, e || [], !1))
        },
        not: function(e) {
          return this.pushStack(M(this, e || [], !0))
        },
        is: function(e) {
          return !!M(this, "string" == typeof e && T.test(e) ? w(e) : e || [], !1).length
        }
      });
      var E, O = /^(?:\s*(<[\w\W]+>)[^>]*|#([\w-]+))$/;
      (w.fn.init = function(e, t, n) {
        var r, i;
        if (!e) return this;
        if (n = n || E, "string" == typeof e) {
          if (!(r = "<" === e[0] && ">" === e[e.length - 1] && e.length >= 3 ? [null, e, null] : O.exec(e)) || !r[1] && t) return !t || t.jquery ? (
            t || n).find(e) : this.constructor(t).find(e);
          if (r[1]) {
            if (t = t instanceof w ? t[0] : t, w.merge(this, w.parseHTML(r[1], t && t.nodeType ? t.ownerDocument || t : m, !0)), D.test(r[1]) && w
              .isPlainObject(t))
              for (r in t) h(this[r]) ? this[r](t[r]) : this.attr(r, t[r]);
            return this
          }
          return (i = m.getElementById(r[2])) && (this[0] = i, this.length = 1), this
        }
        return e.nodeType ? (this[0] = e, this.length = 1, this) : h(e) ? void 0 !== n.ready ? n.ready(e) : e(w) : w.makeArray(e, this)
      }).prototype = w.fn, E = w(m);
      var N = /^(?:parents|prev(?:Until|All))/,
        P = {
          children: !0,
          contents: !0,
          next: !0,
          prev: !0
        };

      function q(e, t) {
        for (;
          (e = e[t]) && 1 !== e.nodeType;);
        return e
      }
      w.fn.extend({
        has: function(e) {
          var t = w(e, this),
            n = t.length;
          return this.filter((function() {
            for (var e = 0; e < n; e++)
              if (w.contains(this, t[e])) return !0
          }))
        },
        closest: function(e, t) {
          var n, r = 0,
            i = this.length,
            o = [],
            a = "string" != typeof e && w(e);
          if (!T.test(e))
            for (; r < i; r++)
              for (n = this[r]; n && n !== t; n = n.parentNode)
                if (n.nodeType < 11 && (a ? a.index(n) > -1 : 1 === n.nodeType && w.find.matchesSelector(n, e))) {
                  o.push(n);
                  break
                } return this.pushStack(o.length > 1 ? w.uniqueSort(o) : o)
        },
        index: function(e) {
          return e ? "string" == typeof e ? s.call(w(e), this[0]) : s.call(this, e.jquery ? e[0] : e) : this[0] && this[0].parentNode ? this
            .first().prevAll().length : -1
        },
        add: function(e, t) {
          return this.pushStack(w.uniqueSort(w.merge(this.get(), w(e, t))))
        },
        addBack: function(e) {
          return this.add(null == e ? this.prevObject : this.prevObject.filter(e))
        }
      }), w.each({
        parent: function(e) {
          var t = e.parentNode;
          return t && 11 !== t.nodeType ? t : null
        },
        parents: function(e) {
          return S(e, "parentNode")
        },
        parentsUntil: function(e, t, n) {
          return S(e, "parentNode", n)
        },
        next: function(e) {
          return q(e, "nextSibling")
        },
        prev: function(e) {
          return q(e, "previousSibling")
        },
        nextAll: function(e) {
          return S(e, "nextSibling")
        },
        prevAll: function(e) {
          return S(e, "previousSibling")
        },
        nextUntil: function(e, t, n) {
          return S(e, "nextSibling", n)
        },
        prevUntil: function(e, t, n) {
          return S(e, "previousSibling", n)
        },
        siblings: function(e) {
          return A((e.parentNode || {}).firstChild, e)
        },
        children: function(e) {
          return A(e.firstChild)
        },
        contents: function(e) {
          return null != e.contentDocument && r(e.contentDocument) ? e.contentDocument : (k(e, "template") && (e = e.content || e), w.merge([],
            e.childNodes))
        }
      }, (function(e, t) {
        w.fn[e] = function(n, r) {
          var i = w.map(this, t, n);
          return "Until" !== e.slice(-5) && (r = n), r && "string" == typeof r && (i = w.filter(r, i)), this.length > 1 && (P[e] || w
            .uniqueSort(i), N.test(e) && i.reverse()), this.pushStack(i)
        }
      }));
      var I = /[^\x20\t\r\n\f]+/g;

      function L(e) {
        return e
      }

      function R(e) {
        throw e
      }

      function j(e, t, n, r) {
        var i;
        try {
          e && h(i = e.promise) ? i.call(e).done(t).fail(n) : e && h(i = e.then) ? i.call(e, t, n) : t.apply(void 0, [e].slice(r))
        } catch (e) {
          n.apply(void 0, [e])
        }
      }
      w.Callbacks = function(e) {
        e = "string" == typeof e ? function(e) {
          var t = {};
          return w.each(e.match(I) || [], (function(e, n) {
            t[n] = !0
          })), t
        }(e) : w.extend({}, e);
        var t, n, r, i, o = [],
          a = [],
          s = -1,
          u = function() {
            for (i = i || e.once, r = t = !0; a.length; s = -1)
              for (n = a.shift(); ++s < o.length;) !1 === o[s].apply(n[0], n[1]) && e.stopOnFalse && (s = o.length, n = !1);
            e.memory || (n = !1), t = !1, i && (o = n ? [] : "")
          },
          l = {
            add: function() {
              return o && (n && !t && (s = o.length - 1, a.push(n)), function t(n) {
                w.each(n, (function(n, r) {
                  h(r) ? e.unique && l.has(r) || o.push(r) : r && r.length && "string" !== b(r) && t(r)
                }))
              }(arguments), n && !t && u()), this
            },
            remove: function() {
              return w.each(arguments, (function(e, t) {
                for (var n;
                  (n = w.inArray(t, o, n)) > -1;) o.splice(n, 1), n <= s && s--
              })), this
            },
            has: function(e) {
              return e ? w.inArray(e, o) > -1 : o.length > 0
            },
            empty: function() {
              return o && (o = []), this
            },
            disable: function() {
              return i = a = [], o = n = "", this
            },
            disabled: function() {
              return !o
            },
            lock: function() {
              return i = a = [], n || t || (o = n = ""), this
            },
            locked: function() {
              return !!i
            },
            fireWith: function(e, n) {
              return i || (n = [e, (n = n || []).slice ? n.slice() : n], a.push(n), t || u()), this
            },
            fire: function() {
              return l.fireWith(this, arguments), this
            },
            fired: function() {
              return !!r
            }
          };
        return l
      }, w.extend({
        Deferred: function(t) {
          var n = [
              ["notify", "progress", w.Callbacks("memory"), w.Callbacks("memory"), 2],
              ["resolve", "done", w.Callbacks("once memory"), w.Callbacks("once memory"), 0, "resolved"],
              ["reject", "fail", w.Callbacks("once memory"), w.Callbacks("once memory"), 1, "rejected"]
            ],
            r = "pending",
            i = {
              state: function() {
                return r
              },
              always: function() {
                return o.done(arguments).fail(arguments), this
              },
              catch: function(e) {
                return i.then(null, e)
              },
              pipe: function() {
                var e = arguments;
                return w.Deferred((function(t) {
                  w.each(n, (function(n, r) {
                    var i = h(e[r[4]]) && e[r[4]];
                    o[r[1]]((function() {
                      var e = i && i.apply(this, arguments);
                      e && h(e.promise) ? e.promise().progress(t.notify).done(t.resolve).fail(t.reject) : t[r[0] + "With"](
                        this, i ? [e] : arguments)
                    }))
                  })), e = null
                })).promise()
              },
              then: function(t, r, i) {
                var o = 0;

                function a(t, n, r, i) {
                  return function() {
                    var s = this,
                      u = arguments,
                      l = function() {
                        var e, l;
                        if (!(t < o)) {
                          if ((e = r.apply(s, u)) === n.promise()) throw new TypeError("Thenable self-resolution");
                          l = e && ("object" == typeof e || "function" == typeof e) && e.then, h(l) ? i ? l.call(e, a(o, n, L, i), a(o, n,
                            R, i)) : (o++, l.call(e, a(o, n, L, i), a(o, n, R, i), a(o, n, L, n.notifyWith))) : (r !== L && (s = void 0,
                            u = [e]), (i || n.resolveWith)(s, u))
                        }
                      },
                      c = i ? l : function() {
                        try {
                          l()
                        } catch (e) {
                          w.Deferred.exceptionHook && w.Deferred.exceptionHook(e, c.stackTrace), t + 1 >= o && (r !== R && (s = void 0,
                            u = [e]), n.rejectWith(s, u))
                        }
                      };
                    t ? c() : (w.Deferred.getStackHook && (c.stackTrace = w.Deferred.getStackHook()), e.setTimeout(c))
                  }
                }
                return w.Deferred((function(e) {
                  n[0][3].add(a(0, e, h(i) ? i : L, e.notifyWith)), n[1][3].add(a(0, e, h(t) ? t : L)), n[2][3].add(a(0, e, h(r) ? r :
                    R))
                })).promise()
              },
              promise: function(e) {
                return null != e ? w.extend(e, i) : i
              }
            },
            o = {};
          return w.each(n, (function(e, t) {
            var a = t[2],
              s = t[5];
            i[t[1]] = a.add, s && a.add((function() {
              r = s
            }), n[3 - e][2].disable, n[3 - e][3].disable, n[0][2].lock, n[0][3].lock), a.add(t[3].fire), o[t[0]] = function() {
              return o[t[0] + "With"](this === o ? void 0 : this, arguments), this
            }, o[t[0] + "With"] = a.fireWith
          })), i.promise(o), t && t.call(o, o), o
        },
        when: function(e) {
          var t = arguments.length,
            n = t,
            r = Array(n),
            o = i.call(arguments),
            a = w.Deferred(),
            s = function(e) {
              return function(n) {
                r[e] = this, o[e] = arguments.length > 1 ? i.call(arguments) : n, --t || a.resolveWith(r, o)
              }
            };
          if (t <= 1 && (j(e, a.done(s(n)).resolve, a.reject, !t), "pending" === a.state() || h(o[n] && o[n].then))) return a.then();
          for (; n--;) j(o[n], s(n), a.reject);
          return a.promise()
        }
      });
      var _ = /^(Eval|Internal|Range|Reference|Syntax|Type|URI)Error$/;
      w.Deferred.exceptionHook = function(t, n) {
        e.console && e.console.warn && t && _.test(t.name) && e.console.warn("jQuery.Deferred exception: " + t.message, t.stack, n)
      }, w.readyException = function(t) {
        e.setTimeout((function() {
          throw t
        }))
      };
      var V = w.Deferred();

      function U() {
        m.removeEventListener("DOMContentLoaded", U), e.removeEventListener("load", U), w.ready()
      }
      w.fn.ready = function(e) {
          return V.then(e).catch((function(e) {
            w.readyException(e)
          })), this
        }, w.extend({
          isReady: !1,
          readyWait: 1,
          ready: function(e) {
            (!0 === e ? --w.readyWait : w.isReady) || (w.isReady = !0, !0 !== e && --w.readyWait > 0 || V.resolveWith(m, [w]))
          }
        }), w.ready.then = V.then, "complete" === m.readyState || "loading" !== m.readyState && !m.documentElement.doScroll ? e.setTimeout(w.ready) :
        (m.addEventListener("DOMContentLoaded", U), e.addEventListener("load", U));
      var H = function(e, t, n, r, i, o, a) {
          var s = 0,
            u = e.length,
            l = null == n;
          if ("object" === b(n))
            for (s in i = !0, n) H(e, t, s, n[s], !0, o, a);
          else if (void 0 !== r && (i = !0, h(r) || (a = !0), l && (a ? (t.call(e, r), t = null) : (l = t, t = function(e, t, n) {
              return l.call(w(e), n)
            })), t))
            for (; s < u; s++) t(e[s], n, a ? r : r.call(e[s], s, t(e[s], n)));
          return i ? e : l ? t.call(e) : u ? t(e[0], n) : o
        },
        F = /^-ms-/,
        B = /-([a-z])/g;

      function G(e, t) {
        return t.toUpperCase()
      }

      function W(e) {
        return e.replace(F, "ms-").replace(B, G)
      }
      var z = function(e) {
        return 1 === e.nodeType || 9 === e.nodeType || !+e.nodeType
      };

      function Y() {
        this.expando = w.expando + Y.uid++
      }
      Y.uid = 1, Y.prototype = {
        cache: function(e) {
          var t = e[this.expando];
          return t || (t = {}, z(e) && (e.nodeType ? e[this.expando] = t : Object.defineProperty(e, this.expando, {
            value: t,
            configurable: !0
          }))), t
        },
        set: function(e, t, n) {
          var r, i = this.cache(e);
          if ("string" == typeof t) i[W(t)] = n;
          else
            for (r in t) i[W(r)] = t[r];
          return i
        },
        get: function(e, t) {
          return void 0 === t ? this.cache(e) : e[this.expando] && e[this.expando][W(t)]
        },
        access: function(e, t, n) {
          return void 0 === t || t && "string" == typeof t && void 0 === n ? this.get(e, t) : (this.set(e, t, n), void 0 !== n ? n : t)
        },
        remove: function(e, t) {
          var n, r = e[this.expando];
          if (void 0 !== r) {
            if (void 0 !== t) {
              n = (t = Array.isArray(t) ? t.map(W) : (t = W(t)) in r ? [t] : t.match(I) || []).length;
              for (; n--;) delete r[t[n]]
            }(void 0 === t || w.isEmptyObject(r)) && (e.nodeType ? e[this.expando] = void 0 : delete e[this.expando])
          }
        },
        hasData: function(e) {
          var t = e[this.expando];
          return void 0 !== t && !w.isEmptyObject(t)
        }
      };
      var K = new Y,
        Q = new Y,
        X = /^(?:\{[\w\W]*\}|\[[\w\W]*\])$/,
        J = /[A-Z]/g;

      function Z(e, t, n) {
        var r;
        if (void 0 === n && 1 === e.nodeType)
          if (r = "data-" + t.replace(J, "-$&").toLowerCase(), "string" == typeof(n = e.getAttribute(r))) {
            try {
              n = function(e) {
                return "true" === e || "false" !== e && ("null" === e ? null : e === +e + "" ? +e : X.test(e) ? JSON.parse(e) : e)
              }(n)
            } catch (e) {}
            Q.set(e, t, n)
          } else n = void 0;
        return n
      }
      w.extend({
        hasData: function(e) {
          return Q.hasData(e) || K.hasData(e)
        },
        data: function(e, t, n) {
          return Q.access(e, t, n)
        },
        removeData: function(e, t) {
          Q.remove(e, t)
        },
        _data: function(e, t, n) {
          return K.access(e, t, n)
        },
        _removeData: function(e, t) {
          K.remove(e, t)
        }
      }), w.fn.extend({
        data: function(e, t) {
          var n, r, i, o = this[0],
            a = o && o.attributes;
          if (void 0 === e) {
            if (this.length && (i = Q.get(o), 1 === o.nodeType && !K.get(o, "hasDataAttrs"))) {
              for (n = a.length; n--;) a[n] && 0 === (r = a[n].name).indexOf("data-") && (r = W(r.slice(5)), Z(o, r, i[r]));
              K.set(o, "hasDataAttrs", !0)
            }
            return i
          }
          return "object" == typeof e ? this.each((function() {
            Q.set(this, e)
          })) : H(this, (function(t) {
            var n;
            if (o && void 0 === t) return void 0 !== (n = Q.get(o, e)) || void 0 !== (n = Z(o, e)) ? n : void 0;
            this.each((function() {
              Q.set(this, e, t)
            }))
          }), null, t, arguments.length > 1, null, !0)
        },
        removeData: function(e) {
          return this.each((function() {
            Q.remove(this, e)
          }))
        }
      }), w.extend({
        queue: function(e, t, n) {
          var r;
          if (e) return t = (t || "fx") + "queue", r = K.get(e, t), n && (!r || Array.isArray(n) ? r = K.access(e, t, w.makeArray(n)) : r.push(
            n)), r || []
        },
        dequeue: function(e, t) {
          t = t || "fx";
          var n = w.queue(e, t),
            r = n.length,
            i = n.shift(),
            o = w._queueHooks(e, t);
          "inprogress" === i && (i = n.shift(), r--), i && ("fx" === t && n.unshift("inprogress"), delete o.stop, i.call(e, (function() {
            w.dequeue(e, t)
          }), o)), !r && o && o.empty.fire()
        },
        _queueHooks: function(e, t) {
          var n = t + "queueHooks";
          return K.get(e, n) || K.access(e, n, {
            empty: w.Callbacks("once memory").add((function() {
              K.remove(e, [t + "queue", n])
            }))
          })
        }
      }), w.fn.extend({
        queue: function(e, t) {
          var n = 2;
          return "string" != typeof e && (t = e, e = "fx", n--), arguments.length < n ? w.queue(this[0], e) : void 0 === t ? this : this.each((
            function() {
              var n = w.queue(this, e, t);
              w._queueHooks(this, e), "fx" === e && "inprogress" !== n[0] && w.dequeue(this, e)
            }))
        },
        dequeue: function(e) {
          return this.each((function() {
            w.dequeue(this, e)
          }))
        },
        clearQueue: function(e) {
          return this.queue(e || "fx", [])
        },
        promise: function(e, t) {
          var n, r = 1,
            i = w.Deferred(),
            o = this,
            a = this.length,
            s = function() {
              --r || i.resolveWith(o, [o])
            };
          for ("string" != typeof e && (t = e, e = void 0), e = e || "fx"; a--;)(n = K.get(o[a], e + "queueHooks")) && n.empty && (r++, n.empty
            .add(s));
          return s(), i.promise(t)
        }
      });
      var ee = /[+-]?(?:\d*\.|)\d+(?:[eE][+-]?\d+|)/.source,
        te = new RegExp("^(?:([+-])=|)(" + ee + ")([a-z%]*)$", "i"),
        ne = ["Top", "Right", "Bottom", "Left"],
        re = m.documentElement,
        ie = function(e) {
          return w.contains(e.ownerDocument, e)
        },
        oe = {
          composed: !0
        };
      re.getRootNode && (ie = function(e) {
        return w.contains(e.ownerDocument, e) || e.getRootNode(oe) === e.ownerDocument
      });
      var ae = function(e, t) {
        return "none" === (e = t || e).style.display || "" === e.style.display && ie(e) && "none" === w.css(e, "display")
      };

      function se(e, t, n, r) {
        var i, o, a = 20,
          s = r ? function() {
            return r.cur()
          } : function() {
            return w.css(e, t, "")
          },
          u = s(),
          l = n && n[3] || (w.cssNumber[t] ? "" : "px"),
          c = e.nodeType && (w.cssNumber[t] || "px" !== l && +u) && te.exec(w.css(e, t));
        if (c && c[3] !== l) {
          for (u /= 2, l = l || c[3], c = +u || 1; a--;) w.style(e, t, c + l), (1 - o) * (1 - (o = s() / u || .5)) <= 0 && (a = 0), c /= o;
          c *= 2, w.style(e, t, c + l), n = n || []
        }
        return n && (c = +c || +u || 0, i = n[1] ? c + (n[1] + 1) * n[2] : +n[2], r && (r.unit = l, r.start = c, r.end = i)), i
      }
      var ue = {};

      function le(e) {
        var t, n = e.ownerDocument,
          r = e.nodeName,
          i = ue[r];
        return i || (t = n.body.appendChild(n.createElement(r)), i = w.css(t, "display"), t.parentNode.removeChild(t), "none" === i && (i = "block"),
          ue[r] = i, i)
      }

      function ce(e, t) {
        for (var n, r, i = [], o = 0, a = e.length; o < a; o++)(r = e[o]).style && (n = r.style.display, t ? ("none" === n && (i[o] = K.get(r,
          "display") || null, i[o] || (r.style.display = "")), "" === r.style.display && ae(r) && (i[o] = le(r))) : "none" !== n && (i[o] =
          "none", K.set(r, "display", n)));
        for (o = 0; o < a; o++) null != i[o] && (e[o].style.display = i[o]);
        return e
      }
      w.fn.extend({
        show: function() {
          return ce(this, !0)
        },
        hide: function() {
          return ce(this)
        },
        toggle: function(e) {
          return "boolean" == typeof e ? e ? this.show() : this.hide() : this.each((function() {
            ae(this) ? w(this).show() : w(this).hide()
          }))
        }
      });
      var de, pe, fe = /^(?:checkbox|radio)$/i,
        he = /<([a-z][^\/\0>\x20\t\r\n\f]*)/i,
        ge = /^$|^module$|\/(?:java|ecma)script/i;
      de = m.createDocumentFragment().appendChild(m.createElement("div")), (pe = m.createElement("input")).setAttribute("type", "radio"), pe
        .setAttribute("checked", "checked"), pe.setAttribute("name", "t"), de.appendChild(pe), f.checkClone = de.cloneNode(!0).cloneNode(!0).lastChild
        .checked, de.innerHTML = "<textarea>x</textarea>", f.noCloneChecked = !!de.cloneNode(!0).lastChild.defaultValue, de.innerHTML =
        "<option></option>", f.option = !!de.lastChild;
      var me = {
        thead: [1, "<table>", "</table>"],
        col: [2, "<table><colgroup>", "</colgroup></table>"],
        tr: [2, "<table><tbody>", "</tbody></table>"],
        td: [3, "<table><tbody><tr>", "</tr></tbody></table>"],
        _default: [0, "", ""]
      };

      function ve(e, t) {
        var n;
        return n = void 0 !== e.getElementsByTagName ? e.getElementsByTagName(t || "*") : void 0 !== e.querySelectorAll ? e.querySelectorAll(t ||
          "*") : [], void 0 === t || t && k(e, t) ? w.merge([e], n) : n
      }

      function $e(e, t) {
        for (var n = 0, r = e.length; n < r; n++) K.set(e[n], "globalEval", !t || K.get(t[n], "globalEval"))
      }
      me.tbody = me.tfoot = me.colgroup = me.caption = me.thead, me.th = me.td, f.option || (me.optgroup = me.option = [1,
        "<select multiple='multiple'>", "</select>"
      ]);
      var be = /<|&#?\w+;/;

      function ye(e, t, n, r, i) {
        for (var o, a, s, u, l, c, d = t.createDocumentFragment(), p = [], f = 0, h = e.length; f < h; f++)
          if ((o = e[f]) || 0 === o)
            if ("object" === b(o)) w.merge(p, o.nodeType ? [o] : o);
            else if (be.test(o)) {
          for (a = a || d.appendChild(t.createElement("div")), s = (he.exec(o) || ["", ""])[1].toLowerCase(), u = me[s] || me._default, a.innerHTML =
            u[1] + w.htmlPrefilter(o) + u[2], c = u[0]; c--;) a = a.lastChild;
          w.merge(p, a.childNodes), (a = d.firstChild).textContent = ""
        } else p.push(t.createTextNode(o));
        for (d.textContent = "", f = 0; o = p[f++];)
          if (r && w.inArray(o, r) > -1) i && i.push(o);
          else if (l = ie(o), a = ve(d.appendChild(o), "script"), l && $e(a), n)
          for (c = 0; o = a[c++];) ge.test(o.type || "") && n.push(o);
        return d
      }
      var we = /^([^.]*)(?:\.(.+)|)/;

      function xe() {
        return !0
      }

      function Ce() {
        return !1
      }

      function Se(e, t) {
        return e === function() {
          try {
            return m.activeElement
          } catch (e) {}
        }() == ("focus" === t)
      }

      function Ae(e, t, n, r, i, o) {
        var a, s;
        if ("object" == typeof t) {
          for (s in "string" != typeof n && (r = r || n, n = void 0), t) Ae(e, s, n, r, t[s], o);
          return e
        }
        if (null == r && null == i ? (i = n, r = n = void 0) : null == i && ("string" == typeof n ? (i = r, r = void 0) : (i = r, r = n, n = void 0)),
          !1 === i) i = Ce;
        else if (!i) return e;
        return 1 === o && (a = i, (i = function(e) {
          return w().off(e), a.apply(this, arguments)
        }).guid = a.guid || (a.guid = w.guid++)), e.each((function() {
          w.event.add(this, t, i, r, n)
        }))
      }

      function Te(e, t, n) {
        n ? (K.set(e, t, !1), w.event.add(e, t, {
          namespace: !1,
          handler: function(e) {
            var r, o, a = K.get(this, t);
            if (1 & e.isTrigger && this[t]) {
              if (a.length)(w.event.special[t] || {}).delegateType && e.stopPropagation();
              else if (a = i.call(arguments), K.set(this, t, a), r = n(this, t), this[t](), a !== (o = K.get(this, t)) || r ? K.set(this, t, !
                  1) : o = {}, a !== o) return e.stopImmediatePropagation(), e.preventDefault(), o && o.value
            } else a.length && (K.set(this, t, {
              value: w.event.trigger(w.extend(a[0], w.Event.prototype), a.slice(1), this)
            }), e.stopImmediatePropagation())
          }
        })) : void 0 === K.get(e, t) && w.event.add(e, t, xe)
      }
      w.event = {
        global: {},
        add: function(e, t, n, r, i) {
          var o, a, s, u, l, c, d, p, f, h, g, m = K.get(e);
          if (z(e))
            for (n.handler && (n = (o = n).handler, i = o.selector), i && w.find.matchesSelector(re, i), n.guid || (n.guid = w.guid++), (u = m
                .events) || (u = m.events = Object.create(null)), (a = m.handle) || (a = m.handle = function(t) {
                return void 0 !== w && w.event.triggered !== t.type ? w.event.dispatch.apply(e, arguments) : void 0
              }), l = (t = (t || "").match(I) || [""]).length; l--;) f = g = (s = we.exec(t[l]) || [])[1], h = (s[2] || "").split(".").sort(),
              f && (d = w.event.special[f] || {}, f = (i ? d.delegateType : d.bindType) || f, d = w.event.special[f] || {}, c = w.extend({
                type: f,
                origType: g,
                data: r,
                handler: n,
                guid: n.guid,
                selector: i,
                needsContext: i && w.expr.match.needsContext.test(i),
                namespace: h.join(".")
              }, o), (p = u[f]) || ((p = u[f] = []).delegateCount = 0, d.setup && !1 !== d.setup.call(e, r, h, a) || e.addEventListener && e
                .addEventListener(f, a)), d.add && (d.add.call(e, c), c.handler.guid || (c.handler.guid = n.guid)), i ? p.splice(p
                .delegateCount++, 0, c) : p.push(c), w.event.global[f] = !0)
        },
        remove: function(e, t, n, r, i) {
          var o, a, s, u, l, c, d, p, f, h, g, m = K.hasData(e) && K.get(e);
          if (m && (u = m.events)) {
            for (l = (t = (t || "").match(I) || [""]).length; l--;)
              if (f = g = (s = we.exec(t[l]) || [])[1], h = (s[2] || "").split(".").sort(), f) {
                for (d = w.event.special[f] || {}, p = u[f = (r ? d.delegateType : d.bindType) || f] || [], s = s[2] && new RegExp("(^|\\.)" + h
                    .join("\\.(?:.*\\.|)") + "(\\.|$)"), a = o = p.length; o--;) c = p[o], !i && g !== c.origType || n && n.guid !== c.guid ||
                  s && !s.test(c.namespace) || r && r !== c.selector && ("**" !== r || !c.selector) || (p.splice(o, 1), c.selector && p
                    .delegateCount--, d.remove && d.remove.call(e, c));
                a && !p.length && (d.teardown && !1 !== d.teardown.call(e, h, m.handle) || w.removeEvent(e, f, m.handle), delete u[f])
              } else
                for (f in u) w.event.remove(e, f + t[l], n, r, !0);
            w.isEmptyObject(u) && K.remove(e, "handle events")
          }
        },
        dispatch: function(e) {
          var t, n, r, i, o, a, s = new Array(arguments.length),
            u = w.event.fix(e),
            l = (K.get(this, "events") || Object.create(null))[u.type] || [],
            c = w.event.special[u.type] || {};
          for (s[0] = u, t = 1; t < arguments.length; t++) s[t] = arguments[t];
          if (u.delegateTarget = this, !c.preDispatch || !1 !== c.preDispatch.call(this, u)) {
            for (a = w.event.handlers.call(this, u, l), t = 0;
              (i = a[t++]) && !u.isPropagationStopped();)
              for (u.currentTarget = i.elem, n = 0;
                (o = i.handlers[n++]) && !u.isImmediatePropagationStopped();) u.rnamespace && !1 !== o.namespace && !u.rnamespace.test(o
                .namespace) || (u.handleObj = o, u.data = o.data, void 0 !== (r = ((w.event.special[o.origType] || {}).handle || o.handler)
                .apply(i.elem, s)) && !1 === (u.result = r) && (u.preventDefault(), u.stopPropagation()));
            return c.postDispatch && c.postDispatch.call(this, u), u.result
          }
        },
        handlers: function(e, t) {
          var n, r, i, o, a, s = [],
            u = t.delegateCount,
            l = e.target;
          if (u && l.nodeType && !("click" === e.type && e.button >= 1))
            for (; l !== this; l = l.parentNode || this)
              if (1 === l.nodeType && ("click" !== e.type || !0 !== l.disabled)) {
                for (o = [], a = {}, n = 0; n < u; n++) void 0 === a[i = (r = t[n]).selector + " "] && (a[i] = r.needsContext ? w(i, this).index(
                  l) > -1 : w.find(i, this, null, [l]).length), a[i] && o.push(r);
                o.length && s.push({
                  elem: l,
                  handlers: o
                })
              } return l = this, u < t.length && s.push({
            elem: l,
            handlers: t.slice(u)
          }), s
        },
        addProp: function(e, t) {
          Object.defineProperty(w.Event.prototype, e, {
            enumerable: !0,
            configurable: !0,
            get: h(t) ? function() {
              if (this.originalEvent) return t(this.originalEvent)
            } : function() {
              if (this.originalEvent) return this.originalEvent[e]
            },
            set: function(t) {
              Object.defineProperty(this, e, {
                enumerable: !0,
                configurable: !0,
                writable: !0,
                value: t
              })
            }
          })
        },
        fix: function(e) {
          return e[w.expando] ? e : new w.Event(e)
        },
        special: {
          load: {
            noBubble: !0
          },
          click: {
            setup: function(e) {
              var t = this || e;
              return fe.test(t.type) && t.click && k(t, "input") && Te(t, "click", xe), !1
            },
            trigger: function(e) {
              var t = this || e;
              return fe.test(t.type) && t.click && k(t, "input") && Te(t, "click"), !0
            },
            _default: function(e) {
              var t = e.target;
              return fe.test(t.type) && t.click && k(t, "input") && K.get(t, "click") || k(t, "a")
            }
          },
          beforeunload: {
            postDispatch: function(e) {
              void 0 !== e.result && e.originalEvent && (e.originalEvent.returnValue = e.result)
            }
          }
        }
      }, w.removeEvent = function(e, t, n) {
        e.removeEventListener && e.removeEventListener(t, n)
      }, w.Event = function(e, t) {
        if (!(this instanceof w.Event)) return new w.Event(e, t);
        e && e.type ? (this.originalEvent = e, this.type = e.type, this.isDefaultPrevented = e.defaultPrevented || void 0 === e.defaultPrevented &&
            !1 === e.returnValue ? xe : Ce, this.target = e.target && 3 === e.target.nodeType ? e.target.parentNode : e.target, this.currentTarget =
            e.currentTarget, this.relatedTarget = e.relatedTarget) : this.type = e, t && w.extend(this, t), this.timeStamp = e && e.timeStamp ||
          Date.now(), this[w.expando] = !0
      }, w.Event.prototype = {
        constructor: w.Event,
        isDefaultPrevented: Ce,
        isPropagationStopped: Ce,
        isImmediatePropagationStopped: Ce,
        isSimulated: !1,
        preventDefault: function() {
          var e = this.originalEvent;
          this.isDefaultPrevented = xe, e && !this.isSimulated && e.preventDefault()
        },
        stopPropagation: function() {
          var e = this.originalEvent;
          this.isPropagationStopped = xe, e && !this.isSimulated && e.stopPropagation()
        },
        stopImmediatePropagation: function() {
          var e = this.originalEvent;
          this.isImmediatePropagationStopped = xe, e && !this.isSimulated && e.stopImmediatePropagation(), this.stopPropagation()
        }
      }, w.each({
        altKey: !0,
        bubbles: !0,
        cancelable: !0,
        changedTouches: !0,
        ctrlKey: !0,
        detail: !0,
        eventPhase: !0,
        metaKey: !0,
        pageX: !0,
        pageY: !0,
        shiftKey: !0,
        view: !0,
        char: !0,
        code: !0,
        charCode: !0,
        key: !0,
        keyCode: !0,
        button: !0,
        buttons: !0,
        clientX: !0,
        clientY: !0,
        offsetX: !0,
        offsetY: !0,
        pointerId: !0,
        pointerType: !0,
        screenX: !0,
        screenY: !0,
        targetTouches: !0,
        toElement: !0,
        touches: !0,
        which: !0
      }, w.event.addProp), w.each({
        focus: "focusin",
        blur: "focusout"
      }, (function(e, t) {
        w.event.special[e] = {
          setup: function() {
            return Te(this, e, Se), !1
          },
          trigger: function() {
            return Te(this, e), !0
          },
          _default: function() {
            return !0
          },
          delegateType: t
        }
      })), w.each({
        mouseenter: "mouseover",
        mouseleave: "mouseout",
        pointerenter: "pointerover",
        pointerleave: "pointerout"
      }, (function(e, t) {
        w.event.special[e] = {
          delegateType: t,
          bindType: t,
          handle: function(e) {
            var n, r = this,
              i = e.relatedTarget,
              o = e.handleObj;
            return i && (i === r || w.contains(r, i)) || (e.type = o.origType, n = o.handler.apply(this, arguments), e.type = t), n
          }
        }
      })), w.fn.extend({
        on: function(e, t, n, r) {
          return Ae(this, e, t, n, r)
        },
        one: function(e, t, n, r) {
          return Ae(this, e, t, n, r, 1)
        },
        off: function(e, t, n) {
          var r, i;
          if (e && e.preventDefault && e.handleObj) return r = e.handleObj, w(e.delegateTarget).off(r.namespace ? r.origType + "." + r
            .namespace : r.origType, r.selector, r.handler), this;
          if ("object" == typeof e) {
            for (i in e) this.off(i, t, e[i]);
            return this
          }
          return !1 !== t && "function" != typeof t || (n = t, t = void 0), !1 === n && (n = Ce), this.each((function() {
            w.event.remove(this, e, n, t)
          }))
        }
      });
      var ke = /<script|<style|<link/i,
        De = /checked\s*(?:[^=]|=\s*.checked.)/i,
        Me = /^\s*<!(?:\[CDATA\[|--)|(?:\]\]|--)>\s*$/g;

      function Ee(e, t) {
        return k(e, "table") && k(11 !== t.nodeType ? t : t.firstChild, "tr") && w(e).children("tbody")[0] || e
      }

      function Oe(e) {
        return e.type = (null !== e.getAttribute("type")) + "/" + e.type, e
      }

      function Ne(e) {
        return "true/" === (e.type || "").slice(0, 5) ? e.type = e.type.slice(5) : e.removeAttribute("type"), e
      }

      function Pe(e, t) {
        var n, r, i, o, a, s;
        if (1 === t.nodeType) {
          if (K.hasData(e) && (s = K.get(e).events))
            for (i in K.remove(t, "handle events"), s)
              for (n = 0, r = s[i].length; n < r; n++) w.event.add(t, i, s[i][n]);
          Q.hasData(e) && (o = Q.access(e), a = w.extend({}, o), Q.set(t, a))
        }
      }

      function qe(e, t) {
        var n = t.nodeName.toLowerCase();
        "input" === n && fe.test(e.type) ? t.checked = e.checked : "input" !== n && "textarea" !== n || (t.defaultValue = e.defaultValue)
      }

      function Ie(e, t, n, r) {
        t = o(t);
        var i, a, s, u, l, c, d = 0,
          p = e.length,
          g = p - 1,
          m = t[0],
          v = h(m);
        if (v || p > 1 && "string" == typeof m && !f.checkClone && De.test(m)) return e.each((function(i) {
          var o = e.eq(i);
          v && (t[0] = m.call(this, i, o.html())), Ie(o, t, n, r)
        }));
        if (p && (a = (i = ye(t, e[0].ownerDocument, !1, e, r)).firstChild, 1 === i.childNodes.length && (i = a), a || r)) {
          for (u = (s = w.map(ve(i, "script"), Oe)).length; d < p; d++) l = i, d !== g && (l = w.clone(l, !0, !0), u && w.merge(s, ve(l, "script"))),
            n.call(e[d], l, d);
          if (u)
            for (c = s[s.length - 1].ownerDocument, w.map(s, Ne), d = 0; d < u; d++) l = s[d], ge.test(l.type || "") && !K.access(l, "globalEval") &&
              w.contains(c, l) && (l.src && "module" !== (l.type || "").toLowerCase() ? w._evalUrl && !l.noModule && w._evalUrl(l.src, {
                nonce: l.nonce || l.getAttribute("nonce")
              }, c) : $(l.textContent.replace(Me, ""), l, c))
        }
        return e
      }

      function Le(e, t, n) {
        for (var r, i = t ? w.filter(t, e) : e, o = 0; null != (r = i[o]); o++) n || 1 !== r.nodeType || w.cleanData(ve(r)), r.parentNode && (n && ie(
          r) && $e(ve(r, "script")), r.parentNode.removeChild(r));
        return e
      }
      w.extend({
        htmlPrefilter: function(e) {
          return e
        },
        clone: function(e, t, n) {
          var r, i, o, a, s = e.cloneNode(!0),
            u = ie(e);
          if (!(f.noCloneChecked || 1 !== e.nodeType && 11 !== e.nodeType || w.isXMLDoc(e)))
            for (a = ve(s), r = 0, i = (o = ve(e)).length; r < i; r++) qe(o[r], a[r]);
          if (t)
            if (n)
              for (o = o || ve(e), a = a || ve(s), r = 0, i = o.length; r < i; r++) Pe(o[r], a[r]);
            else Pe(e, s);
          return (a = ve(s, "script")).length > 0 && $e(a, !u && ve(e, "script")), s
        },
        cleanData: function(e) {
          for (var t, n, r, i = w.event.special, o = 0; void 0 !== (n = e[o]); o++)
            if (z(n)) {
              if (t = n[K.expando]) {
                if (t.events)
                  for (r in t.events) i[r] ? w.event.remove(n, r) : w.removeEvent(n, r, t.handle);
                n[K.expando] = void 0
              }
              n[Q.expando] && (n[Q.expando] = void 0)
            }
        }
      }), w.fn.extend({
        detach: function(e) {
          return Le(this, e, !0)
        },
        remove: function(e) {
          return Le(this, e)
        },
        text: function(e) {
          return H(this, (function(e) {
            return void 0 === e ? w.text(this) : this.empty().each((function() {
              1 !== this.nodeType && 11 !== this.nodeType && 9 !== this.nodeType || (this.textContent = e)
            }))
          }), null, e, arguments.length)
        },
        append: function() {
          return Ie(this, arguments, (function(e) {
            1 !== this.nodeType && 11 !== this.nodeType && 9 !== this.nodeType || Ee(this, e).appendChild(e)
          }))
        },
        prepend: function() {
          return Ie(this, arguments, (function(e) {
            if (1 === this.nodeType || 11 === this.nodeType || 9 === this.nodeType) {
              var t = Ee(this, e);
              t.insertBefore(e, t.firstChild)
            }
          }))
        },
        before: function() {
          return Ie(this, arguments, (function(e) {
            this.parentNode && this.parentNode.insertBefore(e, this)
          }))
        },
        after: function() {
          return Ie(this, arguments, (function(e) {
            this.parentNode && this.parentNode.insertBefore(e, this.nextSibling)
          }))
        },
        empty: function() {
          for (var e, t = 0; null != (e = this[t]); t++) 1 === e.nodeType && (w.cleanData(ve(e, !1)), e.textContent = "");
          return this
        },
        clone: function(e, t) {
          return e = null != e && e, t = null == t ? e : t, this.map((function() {
            return w.clone(this, e, t)
          }))
        },
        html: function(e) {
          return H(this, (function(e) {
            var t = this[0] || {},
              n = 0,
              r = this.length;
            if (void 0 === e && 1 === t.nodeType) return t.innerHTML;
            if ("string" == typeof e && !ke.test(e) && !me[(he.exec(e) || ["", ""])[1].toLowerCase()]) {
              e = w.htmlPrefilter(e);
              try {
                for (; n < r; n++) 1 === (t = this[n] || {}).nodeType && (w.cleanData(ve(t, !1)), t.innerHTML = e);
                t = 0
              } catch (e) {}
            }
            t && this.empty().append(e)
          }), null, e, arguments.length)
        },
        replaceWith: function() {
          var e = [];
          return Ie(this, arguments, (function(t) {
            var n = this.parentNode;
            w.inArray(this, e) < 0 && (w.cleanData(ve(this)), n && n.replaceChild(t, this))
          }), e)
        }
      }), w.each({
        appendTo: "append",
        prependTo: "prepend",
        insertBefore: "before",
        insertAfter: "after",
        replaceAll: "replaceWith"
      }, (function(e, t) {
        w.fn[e] = function(e) {
          for (var n, r = [], i = w(e), o = i.length - 1, s = 0; s <= o; s++) n = s === o ? this : this.clone(!0), w(i[s])[t](n), a.apply(r, n
            .get());
          return this.pushStack(r)
        }
      }));
      var Re = new RegExp("^(" + ee + ")(?!px)[a-z%]+$", "i"),
        je = function(t) {
          var n = t.ownerDocument.defaultView;
          return n && n.opener || (n = e), n.getComputedStyle(t)
        },
        _e = function(e, t, n) {
          var r, i, o = {};
          for (i in t) o[i] = e.style[i], e.style[i] = t[i];
          for (i in r = n.call(e), t) e.style[i] = o[i];
          return r
        },
        Ve = new RegExp(ne.join("|"), "i");

      function Ue(e, t, n) {
        var r, i, o, a, s = e.style;
        return (n = n || je(e)) && ("" !== (a = n.getPropertyValue(t) || n[t]) || ie(e) || (a = w.style(e, t)), !f.pixelBoxStyles() && Re.test(a) &&
          Ve.test(t) && (r = s.width, i = s.minWidth, o = s.maxWidth, s.minWidth = s.maxWidth = s.width = a, a = n.width, s.width = r, s.minWidth =
            i, s.maxWidth = o)), void 0 !== a ? a + "" : a
      }

      function He(e, t) {
        return {
          get: function() {
            if (!e()) return (this.get = t).apply(this, arguments);
            delete this.get
          }
        }
      }! function() {
        function t() {
          if (c) {
            l.style.cssText = "position:absolute;left:-11111px;width:60px;margin-top:1px;padding:0;border:0", c.style.cssText =
              "position:relative;display:block;box-sizing:border-box;overflow:scroll;margin:auto;border:1px;padding:1px;width:60%;top:1%", re
              .appendChild(l).appendChild(c);
            var t = e.getComputedStyle(c);
            r = "1%" !== t.top, u = 12 === n(t.marginLeft), c.style.right = "60%", a = 36 === n(t.right), i = 36 === n(t.width), c.style.position =
              "absolute", o = 12 === n(c.offsetWidth / 3), re.removeChild(l), c = null
          }
        }

        function n(e) {
          return Math.round(parseFloat(e))
        }
        var r, i, o, a, s, u, l = m.createElement("div"),
          c = m.createElement("div");
        c.style && (c.style.backgroundClip = "content-box", c.cloneNode(!0).style.backgroundClip = "", f.clearCloneStyle = "content-box" === c.style
          .backgroundClip, w.extend(f, {
            boxSizingReliable: function() {
              return t(), i
            },
            pixelBoxStyles: function() {
              return t(), a
            },
            pixelPosition: function() {
              return t(), r
            },
            reliableMarginLeft: function() {
              return t(), u
            },
            scrollboxSize: function() {
              return t(), o
            },
            reliableTrDimensions: function() {
              var t, n, r, i;
              return null == s && (t = m.createElement("table"), n = m.createElement("tr"), r = m.createElement("div"), t.style.cssText =
                "position:absolute;left:-11111px;border-collapse:separate", n.style.cssText = "border:1px solid", n.style.height = "1px", r
                .style.height = "9px", r.style.display = "block", re.appendChild(t).appendChild(n).appendChild(r), i = e.getComputedStyle(n),
                s = parseInt(i.height, 10) + parseInt(i.borderTopWidth, 10) + parseInt(i.borderBottomWidth, 10) === n.offsetHeight, re
                .removeChild(t)), s
            }
          }))
      }();
      var Fe = ["Webkit", "Moz", "ms"],
        Be = m.createElement("div").style,
        Ge = {};

      function We(e) {
        var t = w.cssProps[e] || Ge[e];
        return t || (e in Be ? e : Ge[e] = function(e) {
          for (var t = e[0].toUpperCase() + e.slice(1), n = Fe.length; n--;)
            if ((e = Fe[n] + t) in Be) return e
        }(e) || e)
      }
      var ze = /^(none|table(?!-c[ea]).+)/,
        Ye = /^--/,
        Ke = {
          position: "absolute",
          visibility: "hidden",
          display: "block"
        },
        Qe = {
          letterSpacing: "0",
          fontWeight: "400"
        };

      function Xe(e, t, n) {
        var r = te.exec(t);
        return r ? Math.max(0, r[2] - (n || 0)) + (r[3] || "px") : t
      }

      function Je(e, t, n, r, i, o) {
        var a = "width" === t ? 1 : 0,
          s = 0,
          u = 0;
        if (n === (r ? "border" : "content")) return 0;
        for (; a < 4; a += 2) "margin" === n && (u += w.css(e, n + ne[a], !0, i)), r ? ("content" === n && (u -= w.css(e, "padding" + ne[a], !0, i)),
          "margin" !== n && (u -= w.css(e, "border" + ne[a] + "Width", !0, i))) : (u += w.css(e, "padding" + ne[a], !0, i), "padding" !== n ? u += w
          .css(e, "border" + ne[a] + "Width", !0, i) : s += w.css(e, "border" + ne[a] + "Width", !0, i));
        return !r && o >= 0 && (u += Math.max(0, Math.ceil(e["offset" + t[0].toUpperCase() + t.slice(1)] - o - u - s - .5)) || 0), u
      }

      function Ze(e, t, n) {
        var r = je(e),
          i = (!f.boxSizingReliable() || n) && "border-box" === w.css(e, "boxSizing", !1, r),
          o = i,
          a = Ue(e, t, r),
          s = "offset" + t[0].toUpperCase() + t.slice(1);
        if (Re.test(a)) {
          if (!n) return a;
          a = "auto"
        }
        return (!f.boxSizingReliable() && i || !f.reliableTrDimensions() && k(e, "tr") || "auto" === a || !parseFloat(a) && "inline" === w.css(e,
          "display", !1, r)) && e.getClientRects().length && (i = "border-box" === w.css(e, "boxSizing", !1, r), (o = s in e) && (a = e[s])), (a =
          parseFloat(a) || 0) + Je(e, t, n || (i ? "border" : "content"), o, r, a) + "px"
      }

      function et(e, t, n, r, i) {
        return new et.prototype.init(e, t, n, r, i)
      }
      w.extend({
        cssHooks: {
          opacity: {
            get: function(e, t) {
              if (t) {
                var n = Ue(e, "opacity");
                return "" === n ? "1" : n
              }
            }
          }
        },
        cssNumber: {
          animationIterationCount: !0,
          columnCount: !0,
          fillOpacity: !0,
          flexGrow: !0,
          flexShrink: !0,
          fontWeight: !0,
          gridArea: !0,
          gridColumn: !0,
          gridColumnEnd: !0,
          gridColumnStart: !0,
          gridRow: !0,
          gridRowEnd: !0,
          gridRowStart: !0,
          lineHeight: !0,
          opacity: !0,
          order: !0,
          orphans: !0,
          widows: !0,
          zIndex: !0,
          zoom: !0
        },
        cssProps: {},
        style: function(e, t, n, r) {
          if (e && 3 !== e.nodeType && 8 !== e.nodeType && e.style) {
            var i, o, a, s = W(t),
              u = Ye.test(t),
              l = e.style;
            if (u || (t = We(s)), a = w.cssHooks[t] || w.cssHooks[s], void 0 === n) return a && "get" in a && void 0 !== (i = a.get(e, !1, r)) ?
              i : l[t];
            "string" === (o = typeof n) && (i = te.exec(n)) && i[1] && (n = se(e, t, i), o = "number"), null != n && n == n && ("number" !==
              o || u || (n += i && i[3] || (w.cssNumber[s] ? "" : "px")), f.clearCloneStyle || "" !== n || 0 !== t.indexOf("background") || (
                l[t] = "inherit"), a && "set" in a && void 0 === (n = a.set(e, n, r)) || (u ? l.setProperty(t, n) : l[t] = n))
          }
        },
        css: function(e, t, n, r) {
          var i, o, a, s = W(t);
          return Ye.test(t) || (t = We(s)), (a = w.cssHooks[t] || w.cssHooks[s]) && "get" in a && (i = a.get(e, !0, n)), void 0 === i && (i =
              Ue(e, t, r)), "normal" === i && t in Qe && (i = Qe[t]), "" === n || n ? (o = parseFloat(i), !0 === n || isFinite(o) ? o || 0 :
            i) : i
        }
      }), w.each(["height", "width"], (function(e, t) {
        w.cssHooks[t] = {
          get: function(e, n, r) {
            if (n) return !ze.test(w.css(e, "display")) || e.getClientRects().length && e.getBoundingClientRect().width ? Ze(e, t, r) : _e(
              e, Ke, (function() {
                return Ze(e, t, r)
              }))
          },
          set: function(e, n, r) {
            var i, o = je(e),
              a = !f.scrollboxSize() && "absolute" === o.position,
              s = (a || r) && "border-box" === w.css(e, "boxSizing", !1, o),
              u = r ? Je(e, t, r, s, o) : 0;
            return s && a && (u -= Math.ceil(e["offset" + t[0].toUpperCase() + t.slice(1)] - parseFloat(o[t]) - Je(e, t, "border", !1, o) -
              .5)), u && (i = te.exec(n)) && "px" !== (i[3] || "px") && (e.style[t] = n, n = w.css(e, t)), Xe(0, n, u)
          }
        }
      })), w.cssHooks.marginLeft = He(f.reliableMarginLeft, (function(e, t) {
        if (t) return (parseFloat(Ue(e, "marginLeft")) || e.getBoundingClientRect().left - _e(e, {
          marginLeft: 0
        }, (function() {
          return e.getBoundingClientRect().left
        }))) + "px"
      })), w.each({
        margin: "",
        padding: "",
        border: "Width"
      }, (function(e, t) {
        w.cssHooks[e + t] = {
          expand: function(n) {
            for (var r = 0, i = {}, o = "string" == typeof n ? n.split(" ") : [n]; r < 4; r++) i[e + ne[r] + t] = o[r] || o[r - 2] || o[0];
            return i
          }
        }, "margin" !== e && (w.cssHooks[e + t].set = Xe)
      })), w.fn.extend({
        css: function(e, t) {
          return H(this, (function(e, t, n) {
            var r, i, o = {},
              a = 0;
            if (Array.isArray(t)) {
              for (r = je(e), i = t.length; a < i; a++) o[t[a]] = w.css(e, t[a], !1, r);
              return o
            }
            return void 0 !== n ? w.style(e, t, n) : w.css(e, t)
          }), e, t, arguments.length > 1)
        }
      }), w.Tween = et, et.prototype = {
        constructor: et,
        init: function(e, t, n, r, i, o) {
          this.elem = e, this.prop = n, this.easing = i || w.easing._default, this.options = t, this.start = this.now = this.cur(), this.end = r,
            this.unit = o || (w.cssNumber[n] ? "" : "px")
        },
        cur: function() {
          var e = et.propHooks[this.prop];
          return e && e.get ? e.get(this) : et.propHooks._default.get(this)
        },
        run: function(e) {
          var t, n = et.propHooks[this.prop];
          return this.options.duration ? this.pos = t = w.easing[this.easing](e, this.options.duration * e, 0, 1, this.options.duration) : this
            .pos = t = e, this.now = (this.end - this.start) * t + this.start, this.options.step && this.options.step.call(this.elem, this.now,
              this), n && n.set ? n.set(this) : et.propHooks._default.set(this), this
        }
      }, et.prototype.init.prototype = et.prototype, et.propHooks = {
        _default: {
          get: function(e) {
            var t;
            return 1 !== e.elem.nodeType || null != e.elem[e.prop] && null == e.elem.style[e.prop] ? e.elem[e.prop] : (t = w.css(e.elem, e.prop,
              "")) && "auto" !== t ? t : 0
          },
          set: function(e) {
            w.fx.step[e.prop] ? w.fx.step[e.prop](e) : 1 !== e.elem.nodeType || !w.cssHooks[e.prop] && null == e.elem.style[We(e.prop)] ? e.elem[e
              .prop] = e.now : w.style(e.elem, e.prop, e.now + e.unit)
          }
        }
      }, et.propHooks.scrollTop = et.propHooks.scrollLeft = {
        set: function(e) {
          e.elem.nodeType && e.elem.parentNode && (e.elem[e.prop] = e.now)
        }
      }, w.easing = {
        linear: function(e) {
          return e
        },
        swing: function(e) {
          return .5 - Math.cos(e * Math.PI) / 2
        },
        _default: "swing"
      }, w.fx = et.prototype.init, w.fx.step = {};
      var tt, nt, rt = /^(?:toggle|show|hide)$/,
        it = /queueHooks$/;

      function ot() {
        nt && (!1 === m.hidden && e.requestAnimationFrame ? e.requestAnimationFrame(ot) : e.setTimeout(ot, w.fx.interval), w.fx.tick())
      }

      function at() {
        return e.setTimeout((function() {
          tt = void 0
        })), tt = Date.now()
      }

      function st(e, t) {
        var n, r = 0,
          i = {
            height: e
          };
        for (t = t ? 1 : 0; r < 4; r += 2 - t) i["margin" + (n = ne[r])] = i["padding" + n] = e;
        return t && (i.opacity = i.width = e), i
      }

      function ut(e, t, n) {
        for (var r, i = (lt.tweeners[t] || []).concat(lt.tweeners["*"]), o = 0, a = i.length; o < a; o++)
          if (r = i[o].call(n, t, e)) return r
      }

      function lt(e, t, n) {
        var r, i, o = 0,
          a = lt.prefilters.length,
          s = w.Deferred().always((function() {
            delete u.elem
          })),
          u = function() {
            if (i) return !1;
            for (var t = tt || at(), n = Math.max(0, l.startTime + l.duration - t), r = 1 - (n / l.duration || 0), o = 0, a = l.tweens.length; o <
              a; o++) l.tweens[o].run(r);
            return s.notifyWith(e, [l, r, n]), r < 1 && a ? n : (a || s.notifyWith(e, [l, 1, 0]), s.resolveWith(e, [l]), !1)
          },
          l = s.promise({
            elem: e,
            props: w.extend({}, t),
            opts: w.extend(!0, {
              specialEasing: {},
              easing: w.easing._default
            }, n),
            originalProperties: t,
            originalOptions: n,
            startTime: tt || at(),
            duration: n.duration,
            tweens: [],
            createTween: function(t, n) {
              var r = w.Tween(e, l.opts, t, n, l.opts.specialEasing[t] || l.opts.easing);
              return l.tweens.push(r), r
            },
            stop: function(t) {
              var n = 0,
                r = t ? l.tweens.length : 0;
              if (i) return this;
              for (i = !0; n < r; n++) l.tweens[n].run(1);
              return t ? (s.notifyWith(e, [l, 1, 0]), s.resolveWith(e, [l, t])) : s.rejectWith(e, [l, t]), this
            }
          }),
          c = l.props;
        for (! function(e, t) {
            var n, r, i, o, a;
            for (n in e)
              if (i = t[r = W(n)], o = e[n], Array.isArray(o) && (i = o[1], o = e[n] = o[0]), n !== r && (e[r] = o, delete e[n]), (a = w.cssHooks[
                r]) && "expand" in a)
                for (n in o = a.expand(o), delete e[r], o) n in e || (e[n] = o[n], t[n] = i);
              else t[r] = i
          }(c, l.opts.specialEasing); o < a; o++)
          if (r = lt.prefilters[o].call(l, e, c, l.opts)) return h(r.stop) && (w._queueHooks(l.elem, l.opts.queue).stop = r.stop.bind(r)), r;
        return w.map(c, ut, l), h(l.opts.start) && l.opts.start.call(e, l), l.progress(l.opts.progress).done(l.opts.done, l.opts.complete).fail(l.opts
          .fail).always(l.opts.always), w.fx.timer(w.extend(u, {
          elem: e,
          anim: l,
          queue: l.opts.queue
        })), l
      }
      w.Animation = w.extend(lt, {
          tweeners: {
            "*": [function(e, t) {
              var n = this.createTween(e, t);
              return se(n.elem, e, te.exec(t), n), n
            }]
          },
          tweener: function(e, t) {
            h(e) ? (t = e, e = ["*"]) : e = e.match(I);
            for (var n, r = 0, i = e.length; r < i; r++) n = e[r], lt.tweeners[n] = lt.tweeners[n] || [], lt.tweeners[n].unshift(t)
          },
          prefilters: [function(e, t, n) {
            var r, i, o, a, s, u, l, c, d = "width" in t || "height" in t,
              p = this,
              f = {},
              h = e.style,
              g = e.nodeType && ae(e),
              m = K.get(e, "fxshow");
            for (r in n.queue || (null == (a = w._queueHooks(e, "fx")).unqueued && (a.unqueued = 0, s = a.empty.fire, a.empty.fire =
            function() {
                a.unqueued || s()
              }), a.unqueued++, p.always((function() {
                p.always((function() {
                  a.unqueued--, w.queue(e, "fx").length || a.empty.fire()
                }))
              }))), t)
              if (i = t[r], rt.test(i)) {
                if (delete t[r], o = o || "toggle" === i, i === (g ? "hide" : "show")) {
                  if ("show" !== i || !m || void 0 === m[r]) continue;
                  g = !0
                }
                f[r] = m && m[r] || w.style(e, r)
              } if ((u = !w.isEmptyObject(t)) || !w.isEmptyObject(f))
              for (r in d && 1 === e.nodeType && (n.overflow = [h.overflow, h.overflowX, h.overflowY], null == (l = m && m.display) && (l = K
                  .get(e, "display")), "none" === (c = w.css(e, "display")) && (l ? c = l : (ce([e], !0), l = e.style.display || l, c = w.css(
                  e, "display"), ce([e]))), ("inline" === c || "inline-block" === c && null != l) && "none" === w.css(e, "float") && (u || (p
                  .done((function() {
                    h.display = l
                  })), null == l && (c = h.display, l = "none" === c ? "" : c)), h.display = "inline-block")), n.overflow && (h.overflow =
                  "hidden", p.always((function() {
                    h.overflow = n.overflow[0], h.overflowX = n.overflow[1], h.overflowY = n.overflow[2]
                  }))), u = !1, f) u || (m ? "hidden" in m && (g = m.hidden) : m = K.access(e, "fxshow", {
                display: l
              }), o && (m.hidden = !g), g && ce([e], !0), p.done((function() {
                for (r in g || ce([e]), K.remove(e, "fxshow"), f) w.style(e, r, f[r])
              }))), u = ut(g ? m[r] : 0, r, p), r in m || (m[r] = u.start, g && (u.end = u.start, u.start = 0))
          }],
          prefilter: function(e, t) {
            t ? lt.prefilters.unshift(e) : lt.prefilters.push(e)
          }
        }), w.speed = function(e, t, n) {
          var r = e && "object" == typeof e ? w.extend({}, e) : {
            complete: n || !n && t || h(e) && e,
            duration: e,
            easing: n && t || t && !h(t) && t
          };
          return w.fx.off ? r.duration = 0 : "number" != typeof r.duration && (r.duration in w.fx.speeds ? r.duration = w.fx.speeds[r.duration] : r
            .duration = w.fx.speeds._default), null != r.queue && !0 !== r.queue || (r.queue = "fx"), r.old = r.complete, r.complete = function() {
            h(r.old) && r.old.call(this), r.queue && w.dequeue(this, r.queue)
          }, r
        }, w.fn.extend({
          fadeTo: function(e, t, n, r) {
            return this.filter(ae).css("opacity", 0).show().end().animate({
              opacity: t
            }, e, n, r)
          },
          animate: function(e, t, n, r) {
            var i = w.isEmptyObject(e),
              o = w.speed(t, n, r),
              a = function() {
                var t = lt(this, w.extend({}, e), o);
                (i || K.get(this, "finish")) && t.stop(!0)
              };
            return a.finish = a, i || !1 === o.queue ? this.each(a) : this.queue(o.queue, a)
          },
          stop: function(e, t, n) {
            var r = function(e) {
              var t = e.stop;
              delete e.stop, t(n)
            };
            return "string" != typeof e && (n = t, t = e, e = void 0), t && this.queue(e || "fx", []), this.each((function() {
              var t = !0,
                i = null != e && e + "queueHooks",
                o = w.timers,
                a = K.get(this);
              if (i) a[i] && a[i].stop && r(a[i]);
              else
                for (i in a) a[i] && a[i].stop && it.test(i) && r(a[i]);
              for (i = o.length; i--;) o[i].elem !== this || null != e && o[i].queue !== e || (o[i].anim.stop(n), t = !1, o.splice(i, 1));
              !t && n || w.dequeue(this, e)
            }))
          },
          finish: function(e) {
            return !1 !== e && (e = e || "fx"), this.each((function() {
              var t, n = K.get(this),
                r = n[e + "queue"],
                i = n[e + "queueHooks"],
                o = w.timers,
                a = r ? r.length : 0;
              for (n.finish = !0, w.queue(this, e, []), i && i.stop && i.stop.call(this, !0), t = o.length; t--;) o[t].elem === this && o[t]
                .queue === e && (o[t].anim.stop(!0), o.splice(t, 1));
              for (t = 0; t < a; t++) r[t] && r[t].finish && r[t].finish.call(this);
              delete n.finish
            }))
          }
        }), w.each(["toggle", "show", "hide"], (function(e, t) {
          var n = w.fn[t];
          w.fn[t] = function(e, r, i) {
            return null == e || "boolean" == typeof e ? n.apply(this, arguments) : this.animate(st(t, !0), e, r, i)
          }
        })), w.each({
          slideDown: st("show"),
          slideUp: st("hide"),
          slideToggle: st("toggle"),
          fadeIn: {
            opacity: "show"
          },
          fadeOut: {
            opacity: "hide"
          },
          fadeToggle: {
            opacity: "toggle"
          }
        }, (function(e, t) {
          w.fn[e] = function(e, n, r) {
            return this.animate(t, e, n, r)
          }
        })), w.timers = [], w.fx.tick = function() {
          var e, t = 0,
            n = w.timers;
          for (tt = Date.now(); t < n.length; t++)(e = n[t])() || n[t] !== e || n.splice(t--, 1);
          n.length || w.fx.stop(), tt = void 0
        }, w.fx.timer = function(e) {
          w.timers.push(e), w.fx.start()
        }, w.fx.interval = 13, w.fx.start = function() {
          nt || (nt = !0, ot())
        }, w.fx.stop = function() {
          nt = null
        }, w.fx.speeds = {
          slow: 600,
          fast: 200,
          _default: 400
        }, w.fn.delay = function(t, n) {
          return t = w.fx && w.fx.speeds[t] || t, n = n || "fx", this.queue(n, (function(n, r) {
            var i = e.setTimeout(n, t);
            r.stop = function() {
              e.clearTimeout(i)
            }
          }))
        },
        function() {
          var e = m.createElement("input"),
            t = m.createElement("select").appendChild(m.createElement("option"));
          e.type = "checkbox", f.checkOn = "" !== e.value, f.optSelected = t.selected, (e = m.createElement("input")).value = "t", e.type = "radio", f
            .radioValue = "t" === e.value
        }();
      var ct, dt = w.expr.attrHandle;
      w.fn.extend({
        attr: function(e, t) {
          return H(this, w.attr, e, t, arguments.length > 1)
        },
        removeAttr: function(e) {
          return this.each((function() {
            w.removeAttr(this, e)
          }))
        }
      }), w.extend({
        attr: function(e, t, n) {
          var r, i, o = e.nodeType;
          if (3 !== o && 8 !== o && 2 !== o) return void 0 === e.getAttribute ? w.prop(e, t, n) : (1 === o && w.isXMLDoc(e) || (i = w.attrHooks[
              t.toLowerCase()] || (w.expr.match.bool.test(t) ? ct : void 0)), void 0 !== n ? null === n ? void w.removeAttr(e, t) : i &&
            "set" in i && void 0 !== (r = i.set(e, n, t)) ? r : (e.setAttribute(t, n + ""), n) : i && "get" in i && null !== (r = i.get(e,
              t)) ? r : null == (r = w.find.attr(e, t)) ? void 0 : r)
        },
        attrHooks: {
          type: {
            set: function(e, t) {
              if (!f.radioValue && "radio" === t && k(e, "input")) {
                var n = e.value;
                return e.setAttribute("type", t), n && (e.value = n), t
              }
            }
          }
        },
        removeAttr: function(e, t) {
          var n, r = 0,
            i = t && t.match(I);
          if (i && 1 === e.nodeType)
            for (; n = i[r++];) e.removeAttribute(n)
        }
      }), ct = {
        set: function(e, t, n) {
          return !1 === t ? w.removeAttr(e, n) : e.setAttribute(n, n), n
        }
      }, w.each(w.expr.match.bool.source.match(/\w+/g), (function(e, t) {
        var n = dt[t] || w.find.attr;
        dt[t] = function(e, t, r) {
          var i, o, a = t.toLowerCase();
          return r || (o = dt[a], dt[a] = i, i = null != n(e, t, r) ? a : null, dt[a] = o), i
        }
      }));
      var pt = /^(?:input|select|textarea|button)$/i,
        ft = /^(?:a|area)$/i;

      function ht(e) {
        return (e.match(I) || []).join(" ")
      }

      function gt(e) {
        return e.getAttribute && e.getAttribute("class") || ""
      }

      function mt(e) {
        return Array.isArray(e) ? e : "string" == typeof e && e.match(I) || []
      }
      w.fn.extend({
        prop: function(e, t) {
          return H(this, w.prop, e, t, arguments.length > 1)
        },
        removeProp: function(e) {
          return this.each((function() {
            delete this[w.propFix[e] || e]
          }))
        }
      }), w.extend({
        prop: function(e, t, n) {
          var r, i, o = e.nodeType;
          if (3 !== o && 8 !== o && 2 !== o) return 1 === o && w.isXMLDoc(e) || (t = w.propFix[t] || t, i = w.propHooks[t]), void 0 !== n ? i &&
            "set" in i && void 0 !== (r = i.set(e, n, t)) ? r : e[t] = n : i && "get" in i && null !== (r = i.get(e, t)) ? r : e[t]
        },
        propHooks: {
          tabIndex: {
            get: function(e) {
              var t = w.find.attr(e, "tabindex");
              return t ? parseInt(t, 10) : pt.test(e.nodeName) || ft.test(e.nodeName) && e.href ? 0 : -1
            }
          }
        },
        propFix: {
          for: "htmlFor",
          class: "className"
        }
      }), f.optSelected || (w.propHooks.selected = {
        get: function(e) {
          var t = e.parentNode;
          return t && t.parentNode && t.parentNode.selectedIndex, null
        },
        set: function(e) {
          var t = e.parentNode;
          t && (t.selectedIndex, t.parentNode && t.parentNode.selectedIndex)
        }
      }), w.each(["tabIndex", "readOnly", "maxLength", "cellSpacing", "cellPadding", "rowSpan", "colSpan", "useMap", "frameBorder",
        "contentEditable"
      ], (function() {
        w.propFix[this.toLowerCase()] = this
      })), w.fn.extend({
        addClass: function(e) {
          var t, n, r, i, o, a, s, u = 0;
          if (h(e)) return this.each((function(t) {
            w(this).addClass(e.call(this, t, gt(this)))
          }));
          if ((t = mt(e)).length)
            for (; n = this[u++];)
              if (i = gt(n), r = 1 === n.nodeType && " " + ht(i) + " ") {
                for (a = 0; o = t[a++];) r.indexOf(" " + o + " ") < 0 && (r += o + " ");
                i !== (s = ht(r)) && n.setAttribute("class", s)
              } return this
        },
        removeClass: function(e) {
          var t, n, r, i, o, a, s, u = 0;
          if (h(e)) return this.each((function(t) {
            w(this).removeClass(e.call(this, t, gt(this)))
          }));
          if (!arguments.length) return this.attr("class", "");
          if ((t = mt(e)).length)
            for (; n = this[u++];)
              if (i = gt(n), r = 1 === n.nodeType && " " + ht(i) + " ") {
                for (a = 0; o = t[a++];)
                  for (; r.indexOf(" " + o + " ") > -1;) r = r.replace(" " + o + " ", " ");
                i !== (s = ht(r)) && n.setAttribute("class", s)
              } return this
        },
        toggleClass: function(e, t) {
          var n = typeof e,
            r = "string" === n || Array.isArray(e);
          return "boolean" == typeof t && r ? t ? this.addClass(e) : this.removeClass(e) : h(e) ? this.each((function(n) {
            w(this).toggleClass(e.call(this, n, gt(this), t), t)
          })) : this.each((function() {
            var t, i, o, a;
            if (r)
              for (i = 0, o = w(this), a = mt(e); t = a[i++];) o.hasClass(t) ? o.removeClass(t) : o.addClass(t);
            else void 0 !== e && "boolean" !== n || ((t = gt(this)) && K.set(this, "__className__", t), this.setAttribute && this
              .setAttribute("class", t || !1 === e ? "" : K.get(this, "__className__") || ""))
          }))
        },
        hasClass: function(e) {
          var t, n, r = 0;
          for (t = " " + e + " "; n = this[r++];)
            if (1 === n.nodeType && (" " + ht(gt(n)) + " ").indexOf(t) > -1) return !0;
          return !1
        }
      });
      var vt = /\r/g;
      w.fn.extend({
        val: function(e) {
          var t, n, r, i = this[0];
          return arguments.length ? (r = h(e), this.each((function(n) {
              var i;
              1 === this.nodeType && (null == (i = r ? e.call(this, n, w(this).val()) : e) ? i = "" : "number" == typeof i ? i += "" :
                Array.isArray(i) && (i = w.map(i, (function(e) {
                  return null == e ? "" : e + ""
                }))), (t = w.valHooks[this.type] || w.valHooks[this.nodeName.toLowerCase()]) && "set" in t && void 0 !== t.set(this, i,
                  "value") || (this.value = i))
            }))) : i ? (t = w.valHooks[i.type] || w.valHooks[i.nodeName.toLowerCase()]) && "get" in t && void 0 !== (n = t.get(i, "value")) ?
            n : "string" == typeof(n = i.value) ? n.replace(vt, "") : null == n ? "" : n : void 0
        }
      }), w.extend({
        valHooks: {
          option: {
            get: function(e) {
              var t = w.find.attr(e, "value");
              return null != t ? t : ht(w.text(e))
            }
          },
          select: {
            get: function(e) {
              var t, n, r, i = e.options,
                o = e.selectedIndex,
                a = "select-one" === e.type,
                s = a ? null : [],
                u = a ? o + 1 : i.length;
              for (r = o < 0 ? u : a ? o : 0; r < u; r++)
                if (((n = i[r]).selected || r === o) && !n.disabled && (!n.parentNode.disabled || !k(n.parentNode, "optgroup"))) {
                  if (t = w(n).val(), a) return t;
                  s.push(t)
                } return s
            },
            set: function(e, t) {
              for (var n, r, i = e.options, o = w.makeArray(t), a = i.length; a--;)((r = i[a]).selected = w.inArray(w.valHooks.option.get(r),
                o) > -1) && (n = !0);
              return n || (e.selectedIndex = -1), o
            }
          }
        }
      }), w.each(["radio", "checkbox"], (function() {
        w.valHooks[this] = {
          set: function(e, t) {
            if (Array.isArray(t)) return e.checked = w.inArray(w(e).val(), t) > -1
          }
        }, f.checkOn || (w.valHooks[this].get = function(e) {
          return null === e.getAttribute("value") ? "on" : e.value
        })
      })), f.focusin = "onfocusin" in e;
      var $t = /^(?:focusinfocus|focusoutblur)$/,
        bt = function(e) {
          e.stopPropagation()
        };
      w.extend(w.event, {
        trigger: function(t, n, r, i) {
          var o, a, s, u, l, d, p, f, v = [r || m],
            $ = c.call(t, "type") ? t.type : t,
            b = c.call(t, "namespace") ? t.namespace.split(".") : [];
          if (a = f = s = r = r || m, 3 !== r.nodeType && 8 !== r.nodeType && !$t.test($ + w.event.triggered) && ($.indexOf(".") > -1 && (b = $
                .split("."), $ = b.shift(), b.sort()), l = $.indexOf(":") < 0 && "on" + $, (t = t[w.expando] ? t : new w.Event($, "object" ==
                typeof t && t)).isTrigger = i ? 2 : 3, t.namespace = b.join("."), t.rnamespace = t.namespace ? new RegExp("(^|\\.)" + b.join(
                "\\.(?:.*\\.|)") + "(\\.|$)") : null, t.result = void 0, t.target || (t.target = r), n = null == n ? [t] : w.makeArray(n, [t]),
              p = w.event.special[$] || {}, i || !p.trigger || !1 !== p.trigger.apply(r, n))) {
            if (!i && !p.noBubble && !g(r)) {
              for (u = p.delegateType || $, $t.test(u + $) || (a = a.parentNode); a; a = a.parentNode) v.push(a), s = a;
              s === (r.ownerDocument || m) && v.push(s.defaultView || s.parentWindow || e)
            }
            for (o = 0;
              (a = v[o++]) && !t.isPropagationStopped();) f = a, t.type = o > 1 ? u : p.bindType || $, (d = (K.get(a, "events") || Object
              .create(null))[t.type] && K.get(a, "handle")) && d.apply(a, n), (d = l && a[l]) && d.apply && z(a) && (t.result = d.apply(a, n),
              !1 === t.result && t.preventDefault());
            return t.type = $, i || t.isDefaultPrevented() || p._default && !1 !== p._default.apply(v.pop(), n) || !z(r) || l && h(r[$]) && !g(
              r) && ((s = r[l]) && (r[l] = null), w.event.triggered = $, t.isPropagationStopped() && f.addEventListener($, bt), r[$](), t
              .isPropagationStopped() && f.removeEventListener($, bt), w.event.triggered = void 0, s && (r[l] = s)), t.result
          }
        },
        simulate: function(e, t, n) {
          var r = w.extend(new w.Event, n, {
            type: e,
            isSimulated: !0
          });
          w.event.trigger(r, null, t)
        }
      }), w.fn.extend({
        trigger: function(e, t) {
          return this.each((function() {
            w.event.trigger(e, t, this)
          }))
        },
        triggerHandler: function(e, t) {
          var n = this[0];
          if (n) return w.event.trigger(e, t, n, !0)
        }
      }), f.focusin || w.each({
        focus: "focusin",
        blur: "focusout"
      }, (function(e, t) {
        var n = function(e) {
          w.event.simulate(t, e.target, w.event.fix(e))
        };
        w.event.special[t] = {
          setup: function() {
            var r = this.ownerDocument || this.document || this,
              i = K.access(r, t);
            i || r.addEventListener(e, n, !0), K.access(r, t, (i || 0) + 1)
          },
          teardown: function() {
            var r = this.ownerDocument || this.document || this,
              i = K.access(r, t) - 1;
            i ? K.access(r, t, i) : (r.removeEventListener(e, n, !0), K.remove(r, t))
          }
        }
      }));
      var yt = e.location,
        wt = {
          guid: Date.now()
        },
        xt = /\?/;
      w.parseXML = function(t) {
        var n, r;
        if (!t || "string" != typeof t) return null;
        try {
          n = (new e.DOMParser).parseFromString(t, "text/xml")
        } catch (e) {}
        return r = n && n.getElementsByTagName("parsererror")[0], n && !r || w.error("Invalid XML: " + (r ? w.map(r.childNodes, (function(e) {
          return e.textContent
        })).join("\n") : t)), n
      };
      var Ct = /\[\]$/,
        St = /\r?\n/g,
        At = /^(?:submit|button|image|reset|file)$/i,
        Tt = /^(?:input|select|textarea|keygen)/i;

      function kt(e, t, n, r) {
        var i;
        if (Array.isArray(t)) w.each(t, (function(t, i) {
          n || Ct.test(e) ? r(e, i) : kt(e + "[" + ("object" == typeof i && null != i ? t : "") + "]", i, n, r)
        }));
        else if (n || "object" !== b(t)) r(e, t);
        else
          for (i in t) kt(e + "[" + i + "]", t[i], n, r)
      }
      w.param = function(e, t) {
        var n, r = [],
          i = function(e, t) {
            var n = h(t) ? t() : t;
            r[r.length] = encodeURIComponent(e) + "=" + encodeURIComponent(null == n ? "" : n)
          };
        if (null == e) return "";
        if (Array.isArray(e) || e.jquery && !w.isPlainObject(e)) w.each(e, (function() {
          i(this.name, this.value)
        }));
        else
          for (n in e) kt(n, e[n], t, i);
        return r.join("&")
      }, w.fn.extend({
        serialize: function() {
          return w.param(this.serializeArray())
        },
        serializeArray: function() {
          return this.map((function() {
            var e = w.prop(this, "elements");
            return e ? w.makeArray(e) : this
          })).filter((function() {
            var e = this.type;
            return this.name && !w(this).is(":disabled") && Tt.test(this.nodeName) && !At.test(e) && (this.checked || !fe.test(e))
          })).map((function(e, t) {
            var n = w(this).val();
            return null == n ? null : Array.isArray(n) ? w.map(n, (function(e) {
              return {
                name: t.name,
                value: e.replace(St, "\r\n")
              }
            })) : {
              name: t.name,
              value: n.replace(St, "\r\n")
            }
          })).get()
        }
      });
      var Dt = /%20/g,
        Mt = /#.*$/,
        Et = /([?&])_=[^&]*/,
        Ot = /^(.*?):[ \t]*([^\r\n]*)$/gm,
        Nt = /^(?:GET|HEAD)$/,
        Pt = /^\/\//,
        qt = {},
        It = {},
        Lt = "*/".concat("*"),
        Rt = m.createElement("a");

      function jt(e) {
        return function(t, n) {
          "string" != typeof t && (n = t, t = "*");
          var r, i = 0,
            o = t.toLowerCase().match(I) || [];
          if (h(n))
            for (; r = o[i++];) "+" === r[0] ? (r = r.slice(1) || "*", (e[r] = e[r] || []).unshift(n)) : (e[r] = e[r] || []).push(n)
        }
      }

      function _t(e, t, n, r) {
        var i = {},
          o = e === It;

        function a(s) {
          var u;
          return i[s] = !0, w.each(e[s] || [], (function(e, s) {
            var l = s(t, n, r);
            return "string" != typeof l || o || i[l] ? o ? !(u = l) : void 0 : (t.dataTypes.unshift(l), a(l), !1)
          })), u
        }
        return a(t.dataTypes[0]) || !i["*"] && a("*")
      }

      function Vt(e, t) {
        var n, r, i = w.ajaxSettings.flatOptions || {};
        for (n in t) void 0 !== t[n] && ((i[n] ? e : r || (r = {}))[n] = t[n]);
        return r && w.extend(!0, e, r), e
      }
      Rt.href = yt.href, w.extend({
        active: 0,
        lastModified: {},
        etag: {},
        ajaxSettings: {
          url: yt.href,
          type: "GET",
          isLocal: /^(?:about|app|app-storage|.+-extension|file|res|widget):$/.test(yt.protocol),
          global: !0,
          processData: !0,
          async: !0,
          contentType: "application/x-www-form-urlencoded; charset=UTF-8",
          accepts: {
            "*": Lt,
            text: "text/plain",
            html: "text/html",
            xml: "application/xml, text/xml",
            json: "application/json, text/javascript"
          },
          contents: {
            xml: /\bxml\b/,
            html: /\bhtml/,
            json: /\bjson\b/
          },
          responseFields: {
            xml: "responseXML",
            text: "responseText",
            json: "responseJSON"
          },
          converters: {
            "* text": String,
            "text html": !0,
            "text json": JSON.parse,
            "text xml": w.parseXML
          },
          flatOptions: {
            url: !0,
            context: !0
          }
        },
        ajaxSetup: function(e, t) {
          return t ? Vt(Vt(e, w.ajaxSettings), t) : Vt(w.ajaxSettings, e)
        },
        ajaxPrefilter: jt(qt),
        ajaxTransport: jt(It),
        ajax: function(t, n) {
          "object" == typeof t && (n = t, t = void 0), n = n || {};
          var r, i, o, a, s, u, l, c, d, p, f = w.ajaxSetup({}, n),
            h = f.context || f,
            g = f.context && (h.nodeType || h.jquery) ? w(h) : w.event,
            v = w.Deferred(),
            $ = w.Callbacks("once memory"),
            b = f.statusCode || {},
            y = {},
            x = {},
            C = "canceled",
            S = {
              readyState: 0,
              getResponseHeader: function(e) {
                var t;
                if (l) {
                  if (!a)
                    for (a = {}; t = Ot.exec(o);) a[t[1].toLowerCase() + " "] = (a[t[1].toLowerCase() + " "] || []).concat(t[2]);
                  t = a[e.toLowerCase() + " "]
                }
                return null == t ? null : t.join(", ")
              },
              getAllResponseHeaders: function() {
                return l ? o : null
              },
              setRequestHeader: function(e, t) {
                return null == l && (e = x[e.toLowerCase()] = x[e.toLowerCase()] || e, y[e] = t), this
              },
              overrideMimeType: function(e) {
                return null == l && (f.mimeType = e), this
              },
              statusCode: function(e) {
                var t;
                if (e)
                  if (l) S.always(e[S.status]);
                  else
                    for (t in e) b[t] = [b[t], e[t]];
                return this
              },
              abort: function(e) {
                var t = e || C;
                return r && r.abort(t), A(0, t), this
              }
            };
          if (v.promise(S), f.url = ((t || f.url || yt.href) + "").replace(Pt, yt.protocol + "//"), f.type = n.method || n.type || f.method || f
            .type, f.dataTypes = (f.dataType || "*").toLowerCase().match(I) || [""], null == f.crossDomain) {
            u = m.createElement("a");
            try {
              u.href = f.url, u.href = u.href, f.crossDomain = Rt.protocol + "//" + Rt.host != u.protocol + "//" + u.host
            } catch (e) {
              f.crossDomain = !0
            }
          }
          if (f.data && f.processData && "string" != typeof f.data && (f.data = w.param(f.data, f.traditional)), _t(qt, f, n, S), l) return S;
          for (d in (c = w.event && f.global) && 0 == w.active++ && w.event.trigger("ajaxStart"), f.type = f.type.toUpperCase(), f
            .hasContent = !Nt.test(f.type), i = f.url.replace(Mt, ""), f.hasContent ? f.data && f.processData && 0 === (f.contentType || "")
            .indexOf("application/x-www-form-urlencoded") && (f.data = f.data.replace(Dt, "+")) : (p = f.url.slice(i.length), f.data && (f
              .processData || "string" == typeof f.data) && (i += (xt.test(i) ? "&" : "?") + f.data, delete f.data), !1 === f.cache && (i = i
              .replace(Et, "$1"), p = (xt.test(i) ? "&" : "?") + "_=" + wt.guid++ + p), f.url = i + p), f.ifModified && (w.lastModified[i] && S
              .setRequestHeader("If-Modified-Since", w.lastModified[i]), w.etag[i] && S.setRequestHeader("If-None-Match", w.etag[i])), (f
              .data && f.hasContent && !1 !== f.contentType || n.contentType) && S.setRequestHeader("Content-Type", f.contentType), S
            .setRequestHeader("Accept", f.dataTypes[0] && f.accepts[f.dataTypes[0]] ? f.accepts[f.dataTypes[0]] + ("*" !== f.dataTypes[0] ?
              ", " + Lt + "; q=0.01" : "") : f.accepts["*"]), f.headers) S.setRequestHeader(d, f.headers[d]);
          if (f.beforeSend && (!1 === f.beforeSend.call(h, S, f) || l)) return S.abort();
          if (C = "abort", $.add(f.complete), S.done(f.success), S.fail(f.error), r = _t(It, f, n, S)) {
            if (S.readyState = 1, c && g.trigger("ajaxSend", [S, f]), l) return S;
            f.async && f.timeout > 0 && (s = e.setTimeout((function() {
              S.abort("timeout")
            }), f.timeout));
            try {
              l = !1, r.send(y, A)
            } catch (e) {
              if (l) throw e;
              A(-1, e)
            }
          } else A(-1, "No Transport");

          function A(t, n, a, u) {
            var d, p, m, y, x, C = n;
            l || (l = !0, s && e.clearTimeout(s), r = void 0, o = u || "", S.readyState = t > 0 ? 4 : 0, d = t >= 200 && t < 300 || 304 === t,
              a && (y = function(e, t, n) {
                for (var r, i, o, a, s = e.contents, u = e.dataTypes;
                  "*" === u[0];) u.shift(), void 0 === r && (r = e.mimeType || t.getResponseHeader("Content-Type"));
                if (r)
                  for (i in s)
                    if (s[i] && s[i].test(r)) {
                      u.unshift(i);
                      break
                    } if (u[0] in n) o = u[0];
                else {
                  for (i in n) {
                    if (!u[0] || e.converters[i + " " + u[0]]) {
                      o = i;
                      break
                    }
                    a || (a = i)
                  }
                  o = o || a
                }
                if (o) return o !== u[0] && u.unshift(o), n[o]
              }(f, S, a)), !d && w.inArray("script", f.dataTypes) > -1 && w.inArray("json", f.dataTypes) < 0 && (f.converters["text script"] =
                function() {}), y = function(e, t, n, r) {
                var i, o, a, s, u, l = {},
                  c = e.dataTypes.slice();
                if (c[1])
                  for (a in e.converters) l[a.toLowerCase()] = e.converters[a];
                for (o = c.shift(); o;)
                  if (e.responseFields[o] && (n[e.responseFields[o]] = t), !u && r && e.dataFilter && (t = e.dataFilter(t, e.dataType)), u =
                    o, o = c.shift())
                    if ("*" === o) o = u;
                    else if ("*" !== u && u !== o) {
                  if (!(a = l[u + " " + o] || l["* " + o]))
                    for (i in l)
                      if ((s = i.split(" "))[1] === o && (a = l[u + " " + s[0]] || l["* " + s[0]])) {
                        !0 === a ? a = l[i] : !0 !== l[i] && (o = s[0], c.unshift(s[1]));
                        break
                      } if (!0 !== a)
                    if (a && e.throws) t = a(t);
                    else try {
                      t = a(t)
                    } catch (e) {
                      return {
                        state: "parsererror",
                        error: a ? e : "No conversion from " + u + " to " + o
                      }
                    }
                }
                return {
                  state: "success",
                  data: t
                }
              }(f, y, S, d), d ? (f.ifModified && ((x = S.getResponseHeader("Last-Modified")) && (w.lastModified[i] = x), (x = S
                  .getResponseHeader("etag")) && (w.etag[i] = x)), 204 === t || "HEAD" === f.type ? C = "nocontent" : 304 === t ? C =
                "notmodified" : (C = y.state, p = y.data, d = !(m = y.error))) : (m = C, !t && C || (C = "error", t < 0 && (t = 0))), S
              .status = t, S.statusText = (n || C) + "", d ? v.resolveWith(h, [p, C, S]) : v.rejectWith(h, [S, C, m]), S.statusCode(b), b =
              void 0, c && g.trigger(d ? "ajaxSuccess" : "ajaxError", [S, f, d ? p : m]), $.fireWith(h, [S, C]), c && (g.trigger(
                "ajaxComplete", [S, f]), --w.active || w.event.trigger("ajaxStop")))
          }
          return S
        },
        getJSON: function(e, t, n) {
          return w.get(e, t, n, "json")
        },
        getScript: function(e, t) {
          return w.get(e, void 0, t, "script")
        }
      }), w.each(["get", "post"], (function(e, t) {
        w[t] = function(e, n, r, i) {
          return h(n) && (i = i || r, r = n, n = void 0), w.ajax(w.extend({
            url: e,
            type: t,
            dataType: i,
            data: n,
            success: r
          }, w.isPlainObject(e) && e))
        }
      })), w.ajaxPrefilter((function(e) {
        var t;
        for (t in e.headers) "content-type" === t.toLowerCase() && (e.contentType = e.headers[t] || "")
      })), w._evalUrl = function(e, t, n) {
        return w.ajax({
          url: e,
          type: "GET",
          dataType: "script",
          cache: !0,
          async: !1,
          global: !1,
          converters: {
            "text script": function() {}
          },
          dataFilter: function(e) {
            w.globalEval(e, t, n)
          }
        })
      }, w.fn.extend({
        wrapAll: function(e) {
          var t;
          return this[0] && (h(e) && (e = e.call(this[0])), t = w(e, this[0].ownerDocument).eq(0).clone(!0), this[0].parentNode && t
            .insertBefore(this[0]), t.map((function() {
              for (var e = this; e.firstElementChild;) e = e.firstElementChild;
              return e
            })).append(this)), this
        },
        wrapInner: function(e) {
          return h(e) ? this.each((function(t) {
            w(this).wrapInner(e.call(this, t))
          })) : this.each((function() {
            var t = w(this),
              n = t.contents();
            n.length ? n.wrapAll(e) : t.append(e)
          }))
        },
        wrap: function(e) {
          var t = h(e);
          return this.each((function(n) {
            w(this).wrapAll(t ? e.call(this, n) : e)
          }))
        },
        unwrap: function(e) {
          return this.parent(e).not("body").each((function() {
            w(this).replaceWith(this.childNodes)
          })), this
        }
      }), w.expr.pseudos.hidden = function(e) {
        return !w.expr.pseudos.visible(e)
      }, w.expr.pseudos.visible = function(e) {
        return !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length)
      }, w.ajaxSettings.xhr = function() {
        try {
          return new e.XMLHttpRequest
        } catch (e) {}
      };
      var Ut = {
          0: 200,
          1223: 204
        },
        Ht = w.ajaxSettings.xhr();
      f.cors = !!Ht && "withCredentials" in Ht, f.ajax = Ht = !!Ht, w.ajaxTransport((function(t) {
        var n, r;
        if (f.cors || Ht && !t.crossDomain) return {
          send: function(i, o) {
            var a, s = t.xhr();
            if (s.open(t.type, t.url, t.async, t.username, t.password), t.xhrFields)
              for (a in t.xhrFields) s[a] = t.xhrFields[a];
            for (a in t.mimeType && s.overrideMimeType && s.overrideMimeType(t.mimeType), t.crossDomain || i["X-Requested-With"] || (i[
                "X-Requested-With"] = "XMLHttpRequest"), i) s.setRequestHeader(a, i[a]);
            n = function(e) {
                return function() {
                  n && (n = r = s.onload = s.onerror = s.onabort = s.ontimeout = s.onreadystatechange = null, "abort" === e ? s.abort() :
                    "error" === e ? "number" != typeof s.status ? o(0, "error") : o(s.status, s.statusText) : o(Ut[s.status] || s
                      .status, s.statusText, "text" !== (s.responseType || "text") || "string" != typeof s.responseText ? {
                        binary: s.response
                      } : {
                        text: s.responseText
                      }, s.getAllResponseHeaders()))
                }
              }, s.onload = n(), r = s.onerror = s.ontimeout = n("error"), void 0 !== s.onabort ? s.onabort = r : s.onreadystatechange =
              function() {
                4 === s.readyState && e.setTimeout((function() {
                  n && r()
                }))
              }, n = n("abort");
            try {
              s.send(t.hasContent && t.data || null)
            } catch (e) {
              if (n) throw e
            }
          },
          abort: function() {
            n && n()
          }
        }
      })), w.ajaxPrefilter((function(e) {
        e.crossDomain && (e.contents.script = !1)
      })), w.ajaxSetup({
        accepts: {
          script: "text/javascript, application/javascript, application/ecmascript, application/x-ecmascript"
        },
        contents: {
          script: /\b(?:java|ecma)script\b/
        },
        converters: {
          "text script": function(e) {
            return w.globalEval(e), e
          }
        }
      }), w.ajaxPrefilter("script", (function(e) {
        void 0 === e.cache && (e.cache = !1), e.crossDomain && (e.type = "GET")
      })), w.ajaxTransport("script", (function(e) {
        var t, n;
        if (e.crossDomain || e.scriptAttrs) return {
          send: function(r, i) {
            t = w("<script>").attr(e.scriptAttrs || {}).prop({
              charset: e.scriptCharset,
              src: e.url
            }).on("load error", n = function(e) {
              t.remove(), n = null, e && i("error" === e.type ? 404 : 200, e.type)
            }), m.head.appendChild(t[0])
          },
          abort: function() {
            n && n()
          }
        }
      }));
      var Ft, Bt = [],
        Gt = /(=)\?(?=&|$)|\?\?/;
      w.ajaxSetup({
        jsonp: "callback",
        jsonpCallback: function() {
          var e = Bt.pop() || w.expando + "_" + wt.guid++;
          return this[e] = !0, e
        }
      }), w.ajaxPrefilter("json jsonp", (function(t, n, r) {
        var i, o, a, s = !1 !== t.jsonp && (Gt.test(t.url) ? "url" : "string" == typeof t.data && 0 === (t.contentType || "").indexOf(
          "application/x-www-form-urlencoded") && Gt.test(t.data) && "data");
        if (s || "jsonp" === t.dataTypes[0]) return i = t.jsonpCallback = h(t.jsonpCallback) ? t.jsonpCallback() : t.jsonpCallback, s ? t[s] =
          t[s].replace(Gt, "$1" + i) : !1 !== t.jsonp && (t.url += (xt.test(t.url) ? "&" : "?") + t.jsonp + "=" + i), t.converters[
            "script json"] = function() {
            return a || w.error(i + " was not called"), a[0]
          }, t.dataTypes[0] = "json", o = e[i], e[i] = function() {
            a = arguments
          }, r.always((function() {
            void 0 === o ? w(e).removeProp(i) : e[i] = o, t[i] && (t.jsonpCallback = n.jsonpCallback, Bt.push(i)), a && h(o) && o(a[0]),
              a = o = void 0
          })), "script"
      })), f.createHTMLDocument = ((Ft = m.implementation.createHTMLDocument("").body).innerHTML = "<form></form><form></form>", 2 === Ft.childNodes
        .length), w.parseHTML = function(e, t, n) {
        return "string" != typeof e ? [] : ("boolean" == typeof t && (n = t, t = !1), t || (f.createHTMLDocument ? ((r = (t = m.implementation
            .createHTMLDocument("")).createElement("base")).href = m.location.href, t.head.appendChild(r)) : t = m), o = !n && [], (i = D.exec(
          e)) ? [t.createElement(i[1])] : (i = ye([e], t, o), o && o.length && w(o).remove(), w.merge([], i.childNodes)));
        var r, i, o
      }, w.fn.load = function(e, t, n) {
        var r, i, o, a = this,
          s = e.indexOf(" ");
        return s > -1 && (r = ht(e.slice(s)), e = e.slice(0, s)), h(t) ? (n = t, t = void 0) : t && "object" == typeof t && (i = "POST"), a.length >
          0 && w.ajax({
            url: e,
            type: i || "GET",
            dataType: "html",
            data: t
          }).done((function(e) {
            o = arguments, a.html(r ? w("<div>").append(w.parseHTML(e)).find(r) : e)
          })).always(n && function(e, t) {
            a.each((function() {
              n.apply(this, o || [e.responseText, t, e])
            }))
          }), this
      }, w.expr.pseudos.animated = function(e) {
        return w.grep(w.timers, (function(t) {
          return e === t.elem
        })).length
      }, w.offset = {
        setOffset: function(e, t, n) {
          var r, i, o, a, s, u, l = w.css(e, "position"),
            c = w(e),
            d = {};
          "static" === l && (e.style.position = "relative"), s = c.offset(), o = w.css(e, "top"), u = w.css(e, "left"), ("absolute" === l ||
            "fixed" === l) && (o + u).indexOf("auto") > -1 ? (a = (r = c.position()).top, i = r.left) : (a = parseFloat(o) || 0, i = parseFloat(
            u) || 0), h(t) && (t = t.call(e, n, w.extend({}, s))), null != t.top && (d.top = t.top - s.top + a), null != t.left && (d.left = t
            .left - s.left + i), "using" in t ? t.using.call(e, d) : c.css(d)
        }
      }, w.fn.extend({
        offset: function(e) {
          if (arguments.length) return void 0 === e ? this : this.each((function(t) {
            w.offset.setOffset(this, e, t)
          }));
          var t, n, r = this[0];
          return r ? r.getClientRects().length ? (t = r.getBoundingClientRect(), n = r.ownerDocument.defaultView, {
            top: t.top + n.pageYOffset,
            left: t.left + n.pageXOffset
          }) : {
            top: 0,
            left: 0
          } : void 0
        },
        position: function() {
          if (this[0]) {
            var e, t, n, r = this[0],
              i = {
                top: 0,
                left: 0
              };
            if ("fixed" === w.css(r, "position")) t = r.getBoundingClientRect();
            else {
              for (t = this.offset(), n = r.ownerDocument, e = r.offsetParent || n.documentElement; e && (e === n.body || e === n
                  .documentElement) && "static" === w.css(e, "position");) e = e.parentNode;
              e && e !== r && 1 === e.nodeType && ((i = w(e).offset()).top += w.css(e, "borderTopWidth", !0), i.left += w.css(e,
                "borderLeftWidth", !0))
            }
            return {
              top: t.top - i.top - w.css(r, "marginTop", !0),
              left: t.left - i.left - w.css(r, "marginLeft", !0)
            }
          }
        },
        offsetParent: function() {
          return this.map((function() {
            for (var e = this.offsetParent; e && "static" === w.css(e, "position");) e = e.offsetParent;
            return e || re
          }))
        }
      }), w.each({
        scrollLeft: "pageXOffset",
        scrollTop: "pageYOffset"
      }, (function(e, t) {
        var n = "pageYOffset" === t;
        w.fn[e] = function(r) {
          return H(this, (function(e, r, i) {
            var o;
            if (g(e) ? o = e : 9 === e.nodeType && (o = e.defaultView), void 0 === i) return o ? o[t] : e[r];
            o ? o.scrollTo(n ? o.pageXOffset : i, n ? i : o.pageYOffset) : e[r] = i
          }), e, r, arguments.length)
        }
      })), w.each(["top", "left"], (function(e, t) {
        w.cssHooks[t] = He(f.pixelPosition, (function(e, n) {
          if (n) return n = Ue(e, t), Re.test(n) ? w(e).position()[t] + "px" : n
        }))
      })), w.each({
        Height: "height",
        Width: "width"
      }, (function(e, t) {
        w.each({
          padding: "inner" + e,
          content: t,
          "": "outer" + e
        }, (function(n, r) {
          w.fn[r] = function(i, o) {
            var a = arguments.length && (n || "boolean" != typeof i),
              s = n || (!0 === i || !0 === o ? "margin" : "border");
            return H(this, (function(t, n, i) {
              var o;
              return g(t) ? 0 === r.indexOf("outer") ? t["inner" + e] : t.document.documentElement["client" + e] : 9 === t
                .nodeType ? (o = t.documentElement, Math.max(t.body["scroll" + e], o["scroll" + e], t.body["offset" + e], o[
                  "offset" + e], o["client" + e])) : void 0 === i ? w.css(t, n, s) : w.style(t, n, i, s)
            }), t, a ? i : void 0, a)
          }
        }))
      })), w.each(["ajaxStart", "ajaxStop", "ajaxComplete", "ajaxError", "ajaxSuccess", "ajaxSend"], (function(e, t) {
        w.fn[t] = function(e) {
          return this.on(t, e)
        }
      })), w.fn.extend({
        bind: function(e, t, n) {
          return this.on(e, null, t, n)
        },
        unbind: function(e, t) {
          return this.off(e, null, t)
        },
        delegate: function(e, t, n, r) {
          return this.on(t, e, n, r)
        },
        undelegate: function(e, t, n) {
          return 1 === arguments.length ? this.off(e, "**") : this.off(t, e || "**", n)
        },
        hover: function(e, t) {
          return this.mouseenter(e).mouseleave(t || e)
        }
      }), w.each(
        "blur focus focusin focusout resize scroll click dblclick mousedown mouseup mousemove mouseover mouseout mouseenter mouseleave change select submit keydown keypress keyup contextmenu"
        .split(" "), (function(e, t) {
          w.fn[t] = function(e, n) {
            return arguments.length > 0 ? this.on(t, null, e, n) : this.trigger(t)
          }
        }));
      var Wt = /^[\s\uFEFF\xA0]+|[\s\uFEFF\xA0]+$/g;
      w.proxy = function(e, t) {
          var n, r, o;
          if ("string" == typeof t && (n = e[t], t = e, e = n), h(e)) return r = i.call(arguments, 2), (o = function() {
            return e.apply(t || this, r.concat(i.call(arguments)))
          }).guid = e.guid = e.guid || w.guid++, o
        }, w.holdReady = function(e) {
          e ? w.readyWait++ : w.ready(!0)
        }, w.isArray = Array.isArray, w.parseJSON = JSON.parse, w.nodeName = k, w.isFunction = h, w.isWindow = g, w.camelCase = W, w.type = b, w.now =
        Date.now, w.isNumeric = function(e) {
          var t = w.type(e);
          return ("number" === t || "string" === t) && !isNaN(e - parseFloat(e))
        }, w.trim = function(e) {
          return null == e ? "" : (e + "").replace(Wt, "")
        }, "function" == typeof define && define.amd && define("jquery", [], (function() {
          return w
        }));
      var zt = e.jQuery,
        Yt = e.$;
      return w.noConflict = function(t) {
        return e.$ === w && (e.$ = Yt), t && e.jQuery === w && (e.jQuery = zt), w
      }, void 0 === t && (e.jQuery = e.$ = w), w
    }))
  }, {}],
  34: [function(e, t, n) {
    (function(e) {
      (function() {
        ! function(e, r) {
          "object" == typeof n && void 0 !== t ? t.exports = r() : "function" == typeof define && define.amd ? define("underscore", r) : (e =
            "undefined" != typeof globalThis ? globalThis : e || self,
            function() {
              var t = e._,
                n = e._ = r();
              n.noConflict = function() {
                return e._ = t, n
              }
            }())
        }(this, (function() {
          var t = "1.13.1",
            n = "object" == typeof self && self.self === self && self || "object" == typeof e && e.global === e && e || Function("return this")
          () || {},
            r = Array.prototype,
            i = Object.prototype,
            o = "undefined" != typeof Symbol ? Symbol.prototype : null,
            a = r.push,
            s = r.slice,
            u = i.toString,
            l = i.hasOwnProperty,
            c = "undefined" != typeof ArrayBuffer,
            d = "undefined" != typeof DataView,
            p = Array.isArray,
            f = Object.keys,
            h = Object.create,
            g = c && ArrayBuffer.isView,
            m = isNaN,
            v = isFinite,
            $ = !{
              toString: null
            }.propertyIsEnumerable("toString"),
            b = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"],
            y = Math.pow(2, 53) - 1;

          function w(e, t) {
            return t = null == t ? e.length - 1 : +t,
              function() {
                for (var n = Math.max(arguments.length - t, 0), r = Array(n), i = 0; i < n; i++) r[i] = arguments[i + t];
                switch (t) {
                  case 0:
                    return e.call(this, r);
                  case 1:
                    return e.call(this, arguments[0], r);
                  case 2:
                    return e.call(this, arguments[0], arguments[1], r)
                }
                var o = Array(t + 1);
                for (i = 0; i < t; i++) o[i] = arguments[i];
                return o[t] = r, e.apply(this, o)
              }
          }

          function x(e) {
            var t = typeof e;
            return "function" === t || "object" === t && !!e
          }

          function C(e) {
            return void 0 === e
          }

          function S(e) {
            return !0 === e || !1 === e || "[object Boolean]" === u.call(e)
          }

          function A(e) {
            var t = "[object " + e + "]";
            return function(e) {
              return u.call(e) === t
            }
          }
          var T = A("String"),
            k = A("Number"),
            D = A("Date"),
            M = A("RegExp"),
            E = A("Error"),
            O = A("Symbol"),
            N = A("ArrayBuffer"),
            P = A("Function"),
            q = n.document && n.document.childNodes;
          "function" != typeof /./ && "object" != typeof Int8Array && "function" != typeof q && (P = function(e) {
            return "function" == typeof e || !1
          });
          var I = P,
            L = A("Object"),
            R = d && L(new DataView(new ArrayBuffer(8))),
            j = "undefined" != typeof Map && L(new Map),
            _ = A("DataView");
          var V = R ? function(e) {
              return null != e && I(e.getInt8) && N(e.buffer)
            } : _,
            U = p || A("Array");

          function H(e, t) {
            return null != e && l.call(e, t)
          }
          var F = A("Arguments");
          ! function() {
            F(arguments) || (F = function(e) {
              return H(e, "callee")
            })
          }();
          var B = F;

          function G(e) {
            return k(e) && m(e)
          }

          function W(e) {
            return function() {
              return e
            }
          }

          function z(e) {
            return function(t) {
              var n = e(t);
              return "number" == typeof n && n >= 0 && n <= y
            }
          }

          function Y(e) {
            return function(t) {
              return null == t ? void 0 : t[e]
            }
          }
          var K = Y("byteLength"),
            Q = z(K),
            X = /\[object ((I|Ui)nt(8|16|32)|Float(32|64)|Uint8Clamped|Big(I|Ui)nt64)Array\]/;
          var J = c ? function(e) {
              return g ? g(e) && !V(e) : Q(e) && X.test(u.call(e))
            } : W(!1),
            Z = Y("length");

          function ee(e, t) {
            t = function(e) {
              for (var t = {}, n = e.length, r = 0; r < n; ++r) t[e[r]] = !0;
              return {
                contains: function(e) {
                  return t[e]
                },
                push: function(n) {
                  return t[n] = !0, e.push(n)
                }
              }
            }(t);
            var n = b.length,
              r = e.constructor,
              o = I(r) && r.prototype || i,
              a = "constructor";
            for (H(e, a) && !t.contains(a) && t.push(a); n--;)(a = b[n]) in e && e[a] !== o[a] && !t.contains(a) && t.push(a)
          }

          function te(e) {
            if (!x(e)) return [];
            if (f) return f(e);
            var t = [];
            for (var n in e) H(e, n) && t.push(n);
            return $ && ee(e, t), t
          }

          function ne(e, t) {
            var n = te(t),
              r = n.length;
            if (null == e) return !r;
            for (var i = Object(e), o = 0; o < r; o++) {
              var a = n[o];
              if (t[a] !== i[a] || !(a in i)) return !1
            }
            return !0
          }

          function re(e) {
            return e instanceof re ? e : this instanceof re ? void(this._wrapped = e) : new re(e)
          }

          function ie(e) {
            return new Uint8Array(e.buffer || e, e.byteOffset || 0, K(e))
          }
          re.VERSION = t, re.prototype.value = function() {
            return this._wrapped
          }, re.prototype.valueOf = re.prototype.toJSON = re.prototype.value, re.prototype.toString = function() {
            return String(this._wrapped)
          };
          var oe = "[object DataView]";

          function ae(e, t, n, r) {
            if (e === t) return 0 !== e || 1 / e == 1 / t;
            if (null == e || null == t) return !1;
            if (e != e) return t != t;
            var i = typeof e;
            return ("function" === i || "object" === i || "object" == typeof t) && se(e, t, n, r)
          }

          function se(e, t, n, r) {
            e instanceof re && (e = e._wrapped), t instanceof re && (t = t._wrapped);
            var i = u.call(e);
            if (i !== u.call(t)) return !1;
            if (R && "[object Object]" == i && V(e)) {
              if (!V(t)) return !1;
              i = oe
            }
            switch (i) {
              case "[object RegExp]":
              case "[object String]":
                return "" + e == "" + t;
              case "[object Number]":
                return +e != +e ? +t != +t : 0 == +e ? 1 / +e == 1 / t : +e == +t;
              case "[object Date]":
              case "[object Boolean]":
                return +e == +t;
              case "[object Symbol]":
                return o.valueOf.call(e) === o.valueOf.call(t);
              case "[object ArrayBuffer]":
              case oe:
                return se(ie(e), ie(t), n, r)
            }
            var a = "[object Array]" === i;
            if (!a && J(e)) {
              if (K(e) !== K(t)) return !1;
              if (e.buffer === t.buffer && e.byteOffset === t.byteOffset) return !0;
              a = !0
            }
            if (!a) {
              if ("object" != typeof e || "object" != typeof t) return !1;
              var s = e.constructor,
                l = t.constructor;
              if (s !== l && !(I(s) && s instanceof s && I(l) && l instanceof l) && "constructor" in e && "constructor" in t) return !1
            }
            r = r || [];
            for (var c = (n = n || []).length; c--;)
              if (n[c] === e) return r[c] === t;
            if (n.push(e), r.push(t), a) {
              if ((c = e.length) !== t.length) return !1;
              for (; c--;)
                if (!ae(e[c], t[c], n, r)) return !1
            } else {
              var d, p = te(e);
              if (c = p.length, te(t).length !== c) return !1;
              for (; c--;)
                if (!H(t, d = p[c]) || !ae(e[d], t[d], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
          }

          function ue(e) {
            if (!x(e)) return [];
            var t = [];
            for (var n in e) t.push(n);
            return $ && ee(e, t), t
          }

          function le(e) {
            var t = Z(e);
            return function(n) {
              if (null == n) return !1;
              var r = ue(n);
              if (Z(r)) return !1;
              for (var i = 0; i < t; i++)
                if (!I(n[e[i]])) return !1;
              return e !== he || !I(n[ce])
            }
          }
          var ce = "forEach",
            de = ["clear", "delete"],
            pe = ["get", "has", "set"],
            fe = de.concat(ce, pe),
            he = de.concat(pe),
            ge = ["add"].concat(de, ce, "has"),
            me = j ? le(fe) : A("Map"),
            ve = j ? le(he) : A("WeakMap"),
            $e = j ? le(ge) : A("Set"),
            be = A("WeakSet");

          function ye(e) {
            for (var t = te(e), n = t.length, r = Array(n), i = 0; i < n; i++) r[i] = e[t[i]];
            return r
          }

          function we(e) {
            for (var t = {}, n = te(e), r = 0, i = n.length; r < i; r++) t[e[n[r]]] = n[r];
            return t
          }

          function xe(e) {
            var t = [];
            for (var n in e) I(e[n]) && t.push(n);
            return t.sort()
          }

          function Ce(e, t) {
            return function(n) {
              var r = arguments.length;
              if (t && (n = Object(n)), r < 2 || null == n) return n;
              for (var i = 1; i < r; i++)
                for (var o = arguments[i], a = e(o), s = a.length, u = 0; u < s; u++) {
                  var l = a[u];
                  t && void 0 !== n[l] || (n[l] = o[l])
                }
              return n
            }
          }
          var Se = Ce(ue),
            Ae = Ce(te),
            Te = Ce(ue, !0);

          function ke(e) {
            if (!x(e)) return {};
            if (h) return h(e);
            var t = function() {};
            t.prototype = e;
            var n = new t;
            return t.prototype = null, n
          }

          function De(e) {
            return x(e) ? U(e) ? e.slice() : Se({}, e) : e
          }

          function Me(e) {
            return U(e) ? e : [e]
          }

          function Ee(e) {
            return re.toPath(e)
          }

          function Oe(e, t) {
            for (var n = t.length, r = 0; r < n; r++) {
              if (null == e) return;
              e = e[t[r]]
            }
            return n ? e : void 0
          }

          function Ne(e, t, n) {
            var r = Oe(e, Ee(t));
            return C(r) ? n : r
          }

          function Pe(e) {
            return e
          }

          function qe(e) {
            return e = Ae({}, e),
              function(t) {
                return ne(t, e)
              }
          }

          function Ie(e) {
            return e = Ee(e),
              function(t) {
                return Oe(t, e)
              }
          }

          function Le(e, t, n) {
            if (void 0 === t) return e;
            switch (null == n ? 3 : n) {
              case 1:
                return function(n) {
                  return e.call(t, n)
                };
              case 3:
                return function(n, r, i) {
                  return e.call(t, n, r, i)
                };
              case 4:
                return function(n, r, i, o) {
                  return e.call(t, n, r, i, o)
                }
            }
            return function() {
              return e.apply(t, arguments)
            }
          }

          function Re(e, t, n) {
            return null == e ? Pe : I(e) ? Le(e, t, n) : x(e) && !U(e) ? qe(e) : Ie(e)
          }

          function je(e, t) {
            return Re(e, t, 1 / 0)
          }

          function _e(e, t, n) {
            return re.iteratee !== je ? re.iteratee(e, t) : Re(e, t, n)
          }

          function Ve() {}

          function Ue(e, t) {
            return null == t && (t = e, e = 0), e + Math.floor(Math.random() * (t - e + 1))
          }
          re.toPath = Me, re.iteratee = je;
          var He = Date.now || function() {
            return (new Date).getTime()
          };

          function Fe(e) {
            var t = function(t) {
                return e[t]
              },
              n = "(?:" + te(e).join("|") + ")",
              r = RegExp(n),
              i = RegExp(n, "g");
            return function(e) {
              return e = null == e ? "" : "" + e, r.test(e) ? e.replace(i, t) : e
            }
          }
          var Be = {
              "&": "&amp;",
              "<": "&lt;",
              ">": "&gt;",
              '"': "&quot;",
              "'": "&#x27;",
              "`": "&#x60;"
            },
            Ge = Fe(Be),
            We = Fe(we(Be)),
            ze = re.templateSettings = {
              evaluate: /<%([\s\S]+?)%>/g,
              interpolate: /<%=([\s\S]+?)%>/g,
              escape: /<%-([\s\S]+?)%>/g
            },
            Ye = /(.)^/,
            Ke = {
              "'": "'",
              "\\": "\\",
              "\r": "r",
              "\n": "n",
              "\u2028": "u2028",
              "\u2029": "u2029"
            },
            Qe = /\\|'|\r|\n|\u2028|\u2029/g;

          function Xe(e) {
            return "\\" + Ke[e]
          }
          var Je = /^\s*(\w|\$)+\s*$/;
          var Ze = 0;

          function et(e, t, n, r, i) {
            if (!(r instanceof t)) return e.apply(n, i);
            var o = ke(e.prototype),
              a = e.apply(o, i);
            return x(a) ? a : o
          }
          var tt = w((function(e, t) {
            var n = tt.placeholder,
              r = function() {
                for (var i = 0, o = t.length, a = Array(o), s = 0; s < o; s++) a[s] = t[s] === n ? arguments[i++] : t[s];
                for (; i < arguments.length;) a.push(arguments[i++]);
                return et(e, r, this, this, a)
              };
            return r
          }));
          tt.placeholder = re;
          var nt = w((function(e, t, n) {
              if (!I(e)) throw new TypeError("Bind must be called on a function");
              var r = w((function(i) {
                return et(e, r, t, this, n.concat(i))
              }));
              return r
            })),
            rt = z(Z);

          function it(e, t, n, r) {
            if (r = r || [], t || 0 === t) {
              if (t <= 0) return r.concat(e)
            } else t = 1 / 0;
            for (var i = r.length, o = 0, a = Z(e); o < a; o++) {
              var s = e[o];
              if (rt(s) && (U(s) || B(s)))
                if (t > 1) it(s, t - 1, n, r), i = r.length;
                else
                  for (var u = 0, l = s.length; u < l;) r[i++] = s[u++];
              else n || (r[i++] = s)
            }
            return r
          }
          var ot = w((function(e, t) {
            var n = (t = it(t, !1, !1)).length;
            if (n < 1) throw new Error("bindAll must be passed function names");
            for (; n--;) {
              var r = t[n];
              e[r] = nt(e[r], e)
            }
            return e
          }));
          var at = w((function(e, t, n) {
              return setTimeout((function() {
                return e.apply(null, n)
              }), t)
            })),
            st = tt(at, re, 1);

          function ut(e) {
            return function() {
              return !e.apply(this, arguments)
            }
          }

          function lt(e, t) {
            var n;
            return function() {
              return --e > 0 && (n = t.apply(this, arguments)), e <= 1 && (t = null), n
            }
          }
          var ct = tt(lt, 2);

          function dt(e, t, n) {
            t = _e(t, n);
            for (var r, i = te(e), o = 0, a = i.length; o < a; o++)
              if (t(e[r = i[o]], r, e)) return r
          }

          function pt(e) {
            return function(t, n, r) {
              n = _e(n, r);
              for (var i = Z(t), o = e > 0 ? 0 : i - 1; o >= 0 && o < i; o += e)
                if (n(t[o], o, t)) return o;
              return -1
            }
          }
          var ft = pt(1),
            ht = pt(-1);

          function gt(e, t, n, r) {
            for (var i = (n = _e(n, r, 1))(t), o = 0, a = Z(e); o < a;) {
              var s = Math.floor((o + a) / 2);
              n(e[s]) < i ? o = s + 1 : a = s
            }
            return o
          }

          function mt(e, t, n) {
            return function(r, i, o) {
              var a = 0,
                u = Z(r);
              if ("number" == typeof o) e > 0 ? a = o >= 0 ? o : Math.max(o + u, a) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
              else if (n && o && u) return r[o = n(r, i)] === i ? o : -1;
              if (i != i) return (o = t(s.call(r, a, u), G)) >= 0 ? o + a : -1;
              for (o = e > 0 ? a : u - 1; o >= 0 && o < u; o += e)
                if (r[o] === i) return o;
              return -1
            }
          }
          var vt = mt(1, ft, gt),
            $t = mt(-1, ht);

          function bt(e, t, n) {
            var r = (rt(e) ? ft : dt)(e, t, n);
            if (void 0 !== r && -1 !== r) return e[r]
          }

          function yt(e, t, n) {
            var r, i;
            if (t = Le(t, n), rt(e))
              for (r = 0, i = e.length; r < i; r++) t(e[r], r, e);
            else {
              var o = te(e);
              for (r = 0, i = o.length; r < i; r++) t(e[o[r]], o[r], e)
            }
            return e
          }

          function wt(e, t, n) {
            t = _e(t, n);
            for (var r = !rt(e) && te(e), i = (r || e).length, o = Array(i), a = 0; a < i; a++) {
              var s = r ? r[a] : a;
              o[a] = t(e[s], s, e)
            }
            return o
          }

          function xt(e) {
            var t = function(t, n, r, i) {
              var o = !rt(t) && te(t),
                a = (o || t).length,
                s = e > 0 ? 0 : a - 1;
              for (i || (r = t[o ? o[s] : s], s += e); s >= 0 && s < a; s += e) {
                var u = o ? o[s] : s;
                r = n(r, t[u], u, t)
              }
              return r
            };
            return function(e, n, r, i) {
              var o = arguments.length >= 3;
              return t(e, Le(n, i, 4), r, o)
            }
          }
          var Ct = xt(1),
            St = xt(-1);

          function At(e, t, n) {
            var r = [];
            return t = _e(t, n), yt(e, (function(e, n, i) {
              t(e, n, i) && r.push(e)
            })), r
          }

          function Tt(e, t, n) {
            t = _e(t, n);
            for (var r = !rt(e) && te(e), i = (r || e).length, o = 0; o < i; o++) {
              var a = r ? r[o] : o;
              if (!t(e[a], a, e)) return !1
            }
            return !0
          }

          function kt(e, t, n) {
            t = _e(t, n);
            for (var r = !rt(e) && te(e), i = (r || e).length, o = 0; o < i; o++) {
              var a = r ? r[o] : o;
              if (t(e[a], a, e)) return !0
            }
            return !1
          }

          function Dt(e, t, n, r) {
            return rt(e) || (e = ye(e)), ("number" != typeof n || r) && (n = 0), vt(e, t, n) >= 0
          }
          var Mt = w((function(e, t, n) {
            var r, i;
            return I(t) ? i = t : (t = Ee(t), r = t.slice(0, -1), t = t[t.length - 1]), wt(e, (function(e) {
              var o = i;
              if (!o) {
                if (r && r.length && (e = Oe(e, r)), null == e) return;
                o = e[t]
              }
              return null == o ? o : o.apply(e, n)
            }))
          }));

          function Et(e, t) {
            return wt(e, Ie(t))
          }

          function Ot(e, t, n) {
            var r, i, o = -1 / 0,
              a = -1 / 0;
            if (null == t || "number" == typeof t && "object" != typeof e[0] && null != e)
              for (var s = 0, u = (e = rt(e) ? e : ye(e)).length; s < u; s++) null != (r = e[s]) && r > o && (o = r);
            else t = _e(t, n), yt(e, (function(e, n, r) {
              ((i = t(e, n, r)) > a || i === -1 / 0 && o === -1 / 0) && (o = e, a = i)
            }));
            return o
          }

          function Nt(e, t, n) {
            if (null == t || n) return rt(e) || (e = ye(e)), e[Ue(e.length - 1)];
            var r = rt(e) ? De(e) : ye(e),
              i = Z(r);
            t = Math.max(Math.min(t, i), 0);
            for (var o = i - 1, a = 0; a < t; a++) {
              var s = Ue(a, o),
                u = r[a];
              r[a] = r[s], r[s] = u
            }
            return r.slice(0, t)
          }

          function Pt(e, t) {
            return function(n, r, i) {
              var o = t ? [
                [],
                []
              ] : {};
              return r = _e(r, i), yt(n, (function(t, i) {
                var a = r(t, i, n);
                e(o, t, a)
              })), o
            }
          }
          var qt = Pt((function(e, t, n) {
              H(e, n) ? e[n].push(t) : e[n] = [t]
            })),
            It = Pt((function(e, t, n) {
              e[n] = t
            })),
            Lt = Pt((function(e, t, n) {
              H(e, n) ? e[n]++ : e[n] = 1
            })),
            Rt = Pt((function(e, t, n) {
              e[n ? 0 : 1].push(t)
            }), !0),
            jt = /[^\ud800-\udfff]|[\ud800-\udbff][\udc00-\udfff]|[\ud800-\udfff]/g;

          function _t(e, t, n) {
            return t in n
          }
          var Vt = w((function(e, t) {
              var n = {},
                r = t[0];
              if (null == e) return n;
              I(r) ? (t.length > 1 && (r = Le(r, t[1])), t = ue(e)) : (r = _t, t = it(t, !1, !1), e = Object(e));
              for (var i = 0, o = t.length; i < o; i++) {
                var a = t[i],
                  s = e[a];
                r(s, a, e) && (n[a] = s)
              }
              return n
            })),
            Ut = w((function(e, t) {
              var n, r = t[0];
              return I(r) ? (r = ut(r), t.length > 1 && (n = t[1])) : (t = wt(it(t, !1, !1), String), r = function(e, n) {
                return !Dt(t, n)
              }), Vt(e, r, n)
            }));

          function Ht(e, t, n) {
            return s.call(e, 0, Math.max(0, e.length - (null == t || n ? 1 : t)))
          }

          function Ft(e, t, n) {
            return null == e || e.length < 1 ? null == t || n ? void 0 : [] : null == t || n ? e[0] : Ht(e, e.length - t)
          }

          function Bt(e, t, n) {
            return s.call(e, null == t || n ? 1 : t)
          }
          var Gt = w((function(e, t) {
              return t = it(t, !0, !0), At(e, (function(e) {
                return !Dt(t, e)
              }))
            })),
            Wt = w((function(e, t) {
              return Gt(e, t)
            }));

          function zt(e, t, n, r) {
            S(t) || (r = n, n = t, t = !1), null != n && (n = _e(n, r));
            for (var i = [], o = [], a = 0, s = Z(e); a < s; a++) {
              var u = e[a],
                l = n ? n(u, a, e) : u;
              t && !n ? (a && o === l || i.push(u), o = l) : n ? Dt(o, l) || (o.push(l), i.push(u)) : Dt(i, u) || i.push(u)
            }
            return i
          }
          var Yt = w((function(e) {
            return zt(it(e, !0, !0))
          }));

          function Kt(e) {
            for (var t = e && Ot(e, Z).length || 0, n = Array(t), r = 0; r < t; r++) n[r] = Et(e, r);
            return n
          }
          var Qt = w(Kt);

          function Xt(e, t) {
            return e._chain ? re(t).chain() : t
          }

          function Jt(e) {
            return yt(xe(e), (function(t) {
              var n = re[t] = e[t];
              re.prototype[t] = function() {
                var e = [this._wrapped];
                return a.apply(e, arguments), Xt(this, n.apply(re, e))
              }
            })), re
          }
          yt(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], (function(e) {
            var t = r[e];
            re.prototype[e] = function() {
              var n = this._wrapped;
              return null != n && (t.apply(n, arguments), "shift" !== e && "splice" !== e || 0 !== n.length || delete n[0]), Xt(this, n)
            }
          })), yt(["concat", "join", "slice"], (function(e) {
            var t = r[e];
            re.prototype[e] = function() {
              var e = this._wrapped;
              return null != e && (e = t.apply(e, arguments)), Xt(this, e)
            }
          }));
          var Zt = Jt({
            __proto__: null,
            VERSION: t,
            restArguments: w,
            isObject: x,
            isNull: function(e) {
              return null === e
            },
            isUndefined: C,
            isBoolean: S,
            isElement: function(e) {
              return !(!e || 1 !== e.nodeType)
            },
            isString: T,
            isNumber: k,
            isDate: D,
            isRegExp: M,
            isError: E,
            isSymbol: O,
            isArrayBuffer: N,
            isDataView: V,
            isArray: U,
            isFunction: I,
            isArguments: B,
            isFinite: function(e) {
              return !O(e) && v(e) && !isNaN(parseFloat(e))
            },
            isNaN: G,
            isTypedArray: J,
            isEmpty: function(e) {
              if (null == e) return !0;
              var t = Z(e);
              return "number" == typeof t && (U(e) || T(e) || B(e)) ? 0 === t : 0 === Z(te(e))
            },
            isMatch: ne,
            isEqual: function(e, t) {
              return ae(e, t)
            },
            isMap: me,
            isWeakMap: ve,
            isSet: $e,
            isWeakSet: be,
            keys: te,
            allKeys: ue,
            values: ye,
            pairs: function(e) {
              for (var t = te(e), n = t.length, r = Array(n), i = 0; i < n; i++) r[i] = [t[i], e[t[i]]];
              return r
            },
            invert: we,
            functions: xe,
            methods: xe,
            extend: Se,
            extendOwn: Ae,
            assign: Ae,
            defaults: Te,
            create: function(e, t) {
              var n = ke(e);
              return t && Ae(n, t), n
            },
            clone: De,
            tap: function(e, t) {
              return t(e), e
            },
            get: Ne,
            has: function(e, t) {
              for (var n = (t = Ee(t)).length, r = 0; r < n; r++) {
                var i = t[r];
                if (!H(e, i)) return !1;
                e = e[i]
              }
              return !!n
            },
            mapObject: function(e, t, n) {
              t = _e(t, n);
              for (var r = te(e), i = r.length, o = {}, a = 0; a < i; a++) {
                var s = r[a];
                o[s] = t(e[s], s, e)
              }
              return o
            },
            identity: Pe,
            constant: W,
            noop: Ve,
            toPath: Me,
            property: Ie,
            propertyOf: function(e) {
              return null == e ? Ve : function(t) {
                return Ne(e, t)
              }
            },
            matcher: qe,
            matches: qe,
            times: function(e, t, n) {
              var r = Array(Math.max(0, e));
              t = Le(t, n, 1);
              for (var i = 0; i < e; i++) r[i] = t(i);
              return r
            },
            random: Ue,
            now: He,
            escape: Ge,
            unescape: We,
            templateSettings: ze,
            template: function(e, t, n) {
              !t && n && (t = n), t = Te({}, t, re.templateSettings);
              var r = RegExp([(t.escape || Ye).source, (t.interpolate || Ye).source, (t.evaluate || Ye).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
              e.replace(r, (function(t, n, r, a, s) {
                return o += e.slice(i, s).replace(Qe, Xe), i = s + t.length, n ? o += "'+\n((__t=(" + n +
                  "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : a && (o += "';\n" + a +
                    "\n__p+='"), t
              })), o += "';\n";
              var a, s = t.variable;
              if (s) {
                if (!Je.test(s)) throw new Error("variable is not a bare identifier: " + s)
              } else o = "with(obj||{}){\n" + o + "}\n", s = "obj";
              o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
              try {
                a = new Function(s, "_", o)
              } catch (e) {
                throw e.source = o, e
              }
              var u = function(e) {
                return a.call(this, e, re)
              };
              return u.source = "function(" + s + "){\n" + o + "}", u
            },
            result: function(e, t, n) {
              var r = (t = Ee(t)).length;
              if (!r) return I(n) ? n.call(e) : n;
              for (var i = 0; i < r; i++) {
                var o = null == e ? void 0 : e[t[i]];
                void 0 === o && (o = n, i = r), e = I(o) ? o.call(e) : o
              }
              return e
            },
            uniqueId: function(e) {
              var t = ++Ze + "";
              return e ? e + t : t
            },
            chain: function(e) {
              var t = re(e);
              return t._chain = !0, t
            },
            iteratee: je,
            partial: tt,
            bind: nt,
            bindAll: ot,
            memoize: function(e, t) {
              var n = function(r) {
                var i = n.cache,
                  o = "" + (t ? t.apply(this, arguments) : r);
                return H(i, o) || (i[o] = e.apply(this, arguments)), i[o]
              };
              return n.cache = {}, n
            },
            delay: at,
            defer: st,
            throttle: function(e, t, n) {
              var r, i, o, a, s = 0;
              n || (n = {});
              var u = function() {
                  s = !1 === n.leading ? 0 : He(), r = null, a = e.apply(i, o), r || (i = o = null)
                },
                l = function() {
                  var l = He();
                  s || !1 !== n.leading || (s = l);
                  var c = t - (l - s);
                  return i = this, o = arguments, c <= 0 || c > t ? (r && (clearTimeout(r), r = null), s = l, a = e.apply(i, o), r || (i =
                    o = null)) : r || !1 === n.trailing || (r = setTimeout(u, c)), a
                };
              return l.cancel = function() {
                clearTimeout(r), s = 0, r = i = o = null
              }, l
            },
            debounce: function(e, t, n) {
              var r, i, o, a, s, u = function() {
                  var l = He() - i;
                  t > l ? r = setTimeout(u, t - l) : (r = null, n || (a = e.apply(s, o)), r || (o = s = null))
                },
                l = w((function(l) {
                  return s = this, o = l, i = He(), r || (r = setTimeout(u, t), n && (a = e.apply(s, o))), a
                }));
              return l.cancel = function() {
                clearTimeout(r), r = o = s = null
              }, l
            },
            wrap: function(e, t) {
              return tt(t, e)
            },
            negate: ut,
            compose: function() {
              var e = arguments,
                t = e.length - 1;
              return function() {
                for (var n = t, r = e[t].apply(this, arguments); n--;) r = e[n].call(this, r);
                return r
              }
            },
            after: function(e, t) {
              return function() {
                if (--e < 1) return t.apply(this, arguments)
              }
            },
            before: lt,
            once: ct,
            findKey: dt,
            findIndex: ft,
            findLastIndex: ht,
            sortedIndex: gt,
            indexOf: vt,
            lastIndexOf: $t,
            find: bt,
            detect: bt,
            findWhere: function(e, t) {
              return bt(e, qe(t))
            },
            each: yt,
            forEach: yt,
            map: wt,
            collect: wt,
            reduce: Ct,
            foldl: Ct,
            inject: Ct,
            reduceRight: St,
            foldr: St,
            filter: At,
            select: At,
            reject: function(e, t, n) {
              return At(e, ut(_e(t)), n)
            },
            every: Tt,
            all: Tt,
            some: kt,
            any: kt,
            contains: Dt,
            includes: Dt,
            include: Dt,
            invoke: Mt,
            pluck: Et,
            where: function(e, t) {
              return At(e, qe(t))
            },
            max: Ot,
            min: function(e, t, n) {
              var r, i, o = 1 / 0,
                a = 1 / 0;
              if (null == t || "number" == typeof t && "object" != typeof e[0] && null != e)
                for (var s = 0, u = (e = rt(e) ? e : ye(e)).length; s < u; s++) null != (r = e[s]) && r < o && (o = r);
              else t = _e(t, n), yt(e, (function(e, n, r) {
                ((i = t(e, n, r)) < a || i === 1 / 0 && o === 1 / 0) && (o = e, a = i)
              }));
              return o
            },
            shuffle: function(e) {
              return Nt(e, 1 / 0)
            },
            sample: Nt,
            sortBy: function(e, t, n) {
              var r = 0;
              return t = _e(t, n), Et(wt(e, (function(e, n, i) {
                return {
                  value: e,
                  index: r++,
                  criteria: t(e, n, i)
                }
              })).sort((function(e, t) {
                var n = e.criteria,
                  r = t.criteria;
                if (n !== r) {
                  if (n > r || void 0 === n) return 1;
                  if (n < r || void 0 === r) return -1
                }
                return e.index - t.index
              })), "value")
            },
            groupBy: qt,
            indexBy: It,
            countBy: Lt,
            partition: Rt,
            toArray: function(e) {
              return e ? U(e) ? s.call(e) : T(e) ? e.match(jt) : rt(e) ? wt(e, Pe) : ye(e) : []
            },
            size: function(e) {
              return null == e ? 0 : rt(e) ? e.length : te(e).length
            },
            pick: Vt,
            omit: Ut,
            first: Ft,
            head: Ft,
            take: Ft,
            initial: Ht,
            last: function(e, t, n) {
              return null == e || e.length < 1 ? null == t || n ? void 0 : [] : null == t || n ? e[e.length - 1] : Bt(e, Math.max(0, e
                .length - t))
            },
            rest: Bt,
            tail: Bt,
            drop: Bt,
            compact: function(e) {
              return At(e, Boolean)
            },
            flatten: function(e, t) {
              return it(e, t, !1)
            },
            without: Wt,
            uniq: zt,
            unique: zt,
            union: Yt,
            intersection: function(e) {
              for (var t = [], n = arguments.length, r = 0, i = Z(e); r < i; r++) {
                var o = e[r];
                if (!Dt(t, o)) {
                  var a;
                  for (a = 1; a < n && Dt(arguments[a], o); a++);
                  a === n && t.push(o)
                }
              }
              return t
            },
            difference: Gt,
            unzip: Kt,
            transpose: Kt,
            zip: Qt,
            object: function(e, t) {
              for (var n = {}, r = 0, i = Z(e); r < i; r++) t ? n[e[r]] = t[r] : n[e[r][0]] = e[r][1];
              return n
            },
            range: function(e, t, n) {
              null == t && (t = e || 0, e = 0), n || (n = t < e ? -1 : 1);
              for (var r = Math.max(Math.ceil((t - e) / n), 0), i = Array(r), o = 0; o < r; o++, e += n) i[o] = e;
              return i
            },
            chunk: function(e, t) {
              if (null == t || t < 1) return [];
              for (var n = [], r = 0, i = e.length; r < i;) n.push(s.call(e, r, r += t));
              return n
            },
            mixin: Jt,
            default: re
          });
          return Zt._ = Zt, Zt
        }))
      }).call(this)
    }).call(this, "undefined" != typeof global ? global : "undefined" != typeof self ? self : "undefined" != typeof window ? window : {})
  }, {}]
}, {}, [2]);

