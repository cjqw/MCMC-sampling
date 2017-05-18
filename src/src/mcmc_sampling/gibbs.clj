(ns mcmc-sampling.gibbs
  (:require [clojure.core.matrix :as m])
  (:gen-class))

(defn- get-sample
  "TODO: how to get sample?
  MH?"
  [f]
  ([0 0]))

(defn gibbs-get-sample
  "Gibbs algorithm."
  [density-function ]
  (let [partial-x (fn [x] #(density-function [x %]))
        partial-y (fn [y] #(density-function [% y]))
        burn-in 20]
    (loop [i 0
           x 0
           y 0]
      (if (= i burn-in)
        [x y]
        (let [x (get-sample (partial-x x))
              y (get-sample (partial-y y))]
          (recur (inc i) x y))))))
