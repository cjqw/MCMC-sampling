(ns mcmc-sampling.util
  (:gen-class))

(def sqr (fn [x] (* x x)))

(defn uniform-sample
  []
  (rand))

(defn gauss-sample
  "Return a variable which satisfying gauss distribution."
  [mu sigma]
  (let [u (uniform-sample)
        v (uniform-sample)
        x (* (Math/sqrt (* -2 (Math/log u)))
             (Math/cos (* 2 Math/PI v)))]
    (+ mu (* x sigma))))

(defn gauss-function
  [mu sigma x]
  (let [left (/ 1 (Math/sqrt (* 2 Math/PI)))
        index (/ (sqr (- x mu))
                 (* -2 (sqr sigma)))
        right (Math/exp index)]
    (* left right)))
