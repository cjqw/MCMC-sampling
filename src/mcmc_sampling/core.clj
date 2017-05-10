(ns mcmc-sampling.core
  (:require [mcmc-sampling.util :as util]
            [clojure.core.matrix :as m])
  (:gen-class))

(def mu [0 0])
(def a [[1 0] [0 1]])
(def gauss-sample-2d (partial util/gauss-sample-2d mu a))
(def gauss-function-2d (partial util/gauss-function-2d mu a))

(defn- mc-calc-probability-density
  [x d samples]
  (let [neighbor (count (filter (fn [s] (> d (util/distance s x)))
                                samples))
        p (/ neighbor (count samples))
        s (* Math/PI (util/sqr d))]
    (/ p s)))

(defn run
  [x]
  (let [samples (take 1000000 (repeatedly gauss-sample-2d))]
    (mc-calc-probability-density x 0.1 samples)))

(defn -main
  []
  nil)
