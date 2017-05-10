(ns mcmc-sampling.core
  (:require [mcmc-sampling.util :as util]
            [clojure.core.matrix :as m])
  (:gen-class))

(def mu1 [0 0])
(def a1 [[1 0] [0 1]])
(def mu2 [1 1])
(def a2 [[1 2] [3 4]])
(def f1 (partial util/gauss-function-2d mu1 a1))
(def f2 (partial util/gauss-function-2d mu2 a2))
(def standard-gauss-function
  (partial util/gauss-function-2d [0 0] [[1 0] [0 1]]))

(defn density-function
  [x]
  (/ (+ (f1 x) (f2 x)) 2))

(defn get-next-sample
  [previous-sample]
  (let [new-sample (util/gauss-sample-2d)
        threshold (/ (* (standard-gauss-function previous-sample)
                        (density-function new-sample))
                     (* (standard-gauss-function new-sample)
                        (density-function previous-sample)))
        alpha (min threshold 1)
        u (util/uniform-sample)]
    (if (<= u alpha)
      new-sample
      previous-sample)))

(defn- mc-calc-probability-density
  [x d samples]
  (let [neighbor (count (filter (fn [s] (> d (util/distance s x)))
                                samples))
        p (/ neighbor (count samples))
        s (* Math/PI (util/sqr d))]
    (/ p s)))

(defn get-samples
  [n]
  (loop [result [[0 0]]
         i 0]
    (if (= i n)
      result
      (recur (conj result
                   (get-next-sample (peek result)))
             (inc i)))))

(defn get-gauss-samples
  [n]
  (take n (repeatedly util/gauss-sample-2d)))

(defn run
  [x]
  (let [samples (get-samples 100000)]
    (mc-calc-probability-density x 0.1 samples)))

(defn -main
  []
  nil)
