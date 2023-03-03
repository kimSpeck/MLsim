## Roter Faden für die Intro
 - Machine Learning Modelle erfreuen sich großer Beliebtheit, nun auch in der psychologischen Forschung, in der Persönlichkeitspsychologie und insbesondere in der klin. Psychologie
 - Die teilweise in der Forschung berichteten Ergebnisse klingen vielversprechend, Beispiele geben, z.B. Walsh et al. zu suicide prediction, hier auch Literatur zu großen Interaktionseffekten im ML
 - Diese Ergebnisse klingen fast zu gut, um wahr zu sein, und so ist dann auch, Jacobucci et al., 2021
 - Ein Problem scheint dabei die Qualität und Quantität der Daten zu sein, die in der psychologischen Forschung häufig analysiert werden: Geringe sample size (p:n ratio), fuzzy Indikatoren, missingness. All dies behindert die Implementierung von ML-Methoden, die in den Computerwissenschaften entwickelt wurden und die über ein großes N, häufig messfehlerfreie Indikatoren und komplette Datensätze verfügen
- Nach der anfänglicher Euphorie/Hype gibt es zunehmend methodenkritische Betrachtungen, die auf Fehler bei der Umsetzung hinweisen. 1) Kapoor & Narayanan, 2022, information leakage durch falsche Validierung, 2) kein nested resampling oder nur single split, 3) unreliable Indikatoren, Jacobucci & Grimm, 2020, 4) schlechte Praxis wie Validierungssetup oder Hyperparameter-Tuning berichtet werden, Christodoulou et al., 2019
  
## Roter Faden Present Study
 - Mit der Verwendung von ML Algorithmen ist die Hoffnung verbunden, dass man datengetrieben/explorativ Zusammenhänge in den Daten detektiert, ohne diese bereits im Vorfeld spezifiziert zu haben
 - Insbesondere komplexere Modelle wie RF oder GBM sind prinzipiell geeignet nicht-linare und higher-order Effekte zu detektieren
 - ist es erstaunlich, dass nur wenige Studien diese komplexeren Effekte finden, obwohl entsprechende Verfahren verwendet wurden, positives Gegenbeispiel Ali & Ang, 2022, da war die sample size über 3.000 und somit gute Ausgangsbedingungen gegeben
 - Somit stellt sich die Frage, ob typische Datensätze der psychologischen Forschung dafür überhaupt geeignet sind
 - In der vorliegenden Simulationsstudie wollen wir der Frage nachgehen und systematisch die Randbedingungen variieren. Dazu gehören ... s. 2.
 
 ## Forschungsfragen
 1. Unter welchen Randbedingungen sind nicht-lineare bzw. higher-order effects überhaupt detektierbar? [Offene Frage: Wie quantifiziert man den Grad des Auffindens]
 2. Unter welchen Randbedingungen zeigt sich ein Vorteil von GBM gegenüber regularisierten, linearen Modellen.
 3. Lassen sich diese potenziellen Vorteile eines GBM gegenüber linreg auch dann noch finden, wenn die higher-order interaction effects vorab spezifiziert wurden.  [Offene Frage: Bis zu welcher Tiefe?]
 
## 2. Welche Faktoren sollen manipuliert werden? Welche Stufen?
 - sample size, vielleicht einfach p fix halten und so das n:p ratio variieren [100, 200, 500, 1000, 2000]
 - Form der Interaktion, Stärke der Interaktion
 - Reliabilität der Indikatoren [.6, .8, 1]
 - Missingness auf den Indikatoren (?) [5%, 10%, 20%]
