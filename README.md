# MLsim

## Protokoll vom 02.02.:

- Wie wählt man die Korrelationen der Prädiktoren?
- Wie wählt man die Regressionskoeffizienten in der Population?

- generell weniger Prädiktoren -> weg von Big Data?

```Life Satisfaction = age + age^2 + mental_health + physical_health + age*physical_health + age*social_support + mental_health*social_support```
https://psyarxiv.com/gupxj/download

- 3 Bedingungen: alle Macht den Interaktionen, alle Macht den einfachen Effekten, alles gleichmäßig verteilen
- erstmal alle Korrelationen auf 0.3 setzen
(**aktueller Stand**)
  - ~~KS: GitHub-Projekt anlegen~~
  - ~~KJ: lädt den Code hoch~~
  - ~~FS: simpleres Modell simulieren~~
  - **KJ: Modelle testen**
  - **KS: Simulationscode checken**
  - US: Recherche zu den wahren Effekten
  - KS: wie programmiert man die Prädiktorkorrelationsmatrix sinnvoll

17.02.23:
- Was ist mit kategorialen Variablen?
- Welche Rolle spielt es, was man dem enet anbietet? (z.B. mit vs. ohne Interaktionen)
- Wann haben die Bäume echte, genuine Vorteile?
- Was ist mit anderen ML-Modellen?

# Homework
## 1. Forschungsfragen verschriftlichen (US/KJ)
 - Unter welchen Randbedingungen sind nicht-lineare bzw. higher-order effects überhaupt detektierbar?
 - Wann zeigt sich ein Vorteil von GBM gegenüber regularisierten, linearen Modellen bzw. einem mit higher-order interaction.
 [Historie wo kommt man hier, was will man zeigen?]
 
## 2. Welche Faktoren sollen manipuliert werden? Welche Stufen? (alle)
 - sample size (100, 300, 1.000, 10.000)
 - Form der Interaktion (doppelte Zweifachinteraktion, quadratischer Zusammenhang)
 - Höhe der Effekte
 - Aufteilung der Effekte (linear vs. Interaktion)
 - R²
 - ...
 
## 3. Datensatz erstellen mit zwei unterschiedlichen grundlegenden Modellen (KJ)
  - Age > 0 vs. Age <= 0
  
## 4. Algorithmus anpassen, so dass binäre Daten produziert werden (KS)

## 5. R^2 Regressionskoeffizienten zwischen den Extremen (FS)
