###############################################################################
# Choice-test en cages : préférence de Lygus entre le concombre et quatre
# plantes alternatives (luzerne, moutarde, sarrasin, tournesol) ou un témoin
# négatif (pot de terreau sans plante).
#
# Analyse : modèle linéaire mixte (ANOVA à deux facteurs fixes, traitement
#           et position du mur, et un facteur aléatoire, le jour) sur la
#           proportion de Lygus trouvées dans la cage du concombre ;
#           moyennes ajustées comparées à 0.5 ; comparaisons multiples
#           (Tukey, Dunnett) ; analyse complémentaire de la configuration
#           des cages dans le traitement tournesol ; figures.
#
# Données : resultats_cages.csv (séparateur virgule), dans le dossier de
#           travail. Une ligne par essai ; colonnes :
#           id, date, ttt (traitement : N, L, M, S, T),
#           mur (N = aucun ; A = côté plante alternative ; C = côté concombre),
#           total (Lygus lâchées), C (cage concombre), A (cage plante
#           alternative), T (tunnel), M (mortes), I (introuvables),
#           grd_cage (traitement tournesol uniquement : A = seul le
#           tournesol en grande cage ; AC = tournesol et concombre en
#           grandes cages ; N/A pour les autres traitements).
#
# Packages : lme4, lmerTest, emmeans, multcomp, multcompView, ggplot2,
#            patchwork.
###############################################################################

library(lmerTest)   # charge aussi lme4
library(emmeans)
library(multcomp)
library(ggplot2)
library(patchwork)  # assemblage des panneaux de la figure 2


# 1. DONNÉES ##################################################################

d <- read.csv("resultats_cages.csv", fileEncoding = "UTF-8",
              stringsAsFactors = FALSE)

# "T" (tunnel) est renommé pour éviter toute confusion avec TRUE.
names(d)[names(d) == "T"] <- "Tunnel"
names(d)[names(d) == "M"] <- "Mortes"
names(d)[names(d) == "I"] <- "Introuvables"

d$ttt <- factor(d$ttt, levels = c("N", "L", "M", "S", "T"))
d$mur <- factor(d$mur, levels = c("N", "A", "C"))
d$jour <- factor(d$date)

# Variable analysée : proportion de Lygus dans la cage du concombre parmi
# celles ayant choisi une cage. Les Lygus restées dans le tunnel, mortes ou
# introuvables n'ont pas exprimé de choix et sont exclues de ce calcul.
d$choix        <- d$C + d$A
d$prop_C       <- d$C / d$choix
d$taux_reponse <- d$choix / d$total

stopifnot(all(d$C + d$A + d$Tunnel + d$Mortes + d$Introuvables == d$total))

# Les trois essais pilotes (5 Lygus au lieu de 10) sont exclus.
dm <- droplevels(subset(d, total == 10))

# Plan d'expérience : essais par traitement et par position du mur
addmargins(table(Traitement = dm$ttt, Mur = dm$mur))

# Proportions brutes par traitement (non corrigées de l'effet du mur)
aggregate(prop_C ~ ttt, data = dm,
          FUN = function(x) c(n = length(x), moyenne = mean(x), et = sd(x)))


# 2. MODÈLE LINÉAIRE MIXTE ###################################################
# Facteurs fixes : traitement et mur. Le mur (variable parasite) est inclus
# pour corriger son effet plutôt que d'écarter les essais concernés.
# Facteur aléatoire : jour. Les essais d'un même jour partagent les mêmes
# conditions, et le témoin n'a été testé que la dernière semaine.
# Modèle additif, sans interaction (hypothèse vérifiée en 3.3). Le plan
# n'étant pas parfaitement équilibré, chaque facteur fixe est testé après
# ajustement pour l'autre (sommes des carrés de type II), avec les degrés de
# liberté de Satterthwaite.

m <- lmer(prop_C ~ ttt + mur + (1 | jour), data = dm)
anova(m, type = 2)

# Effet aléatoire du jour : variance estimée et test du rapport de
# vraisemblance. Une variance estimée à 0 s'accompagne du message
# "boundary (singular) fit" : les jours n'ajoutent alors aucune variation
# et le modèle se confond avec l'ANOVA à deux facteurs sans le jour.
VarCorr(m)
ranova(m)

# Modèle équivalent sans le jour (identique au modèle mixte lorsque la
# variance du jour est nulle) : sert au R² et à la distance de Cook (3.5).
a <- lm(prop_C ~ ttt + mur, data = dm)
summary(a)$r.squared


# 3. CONDITIONS D'APPLICATION #################################################
# L'indépendance des observations découle du protocole (chaque Lygus utilisée
# une seule fois, plantes tirées au sort, un essai = une unité
# expérimentale).

# 3.1 Normalité des résidus
shapiro.test(residuals(m))

# 3.2 Homogénéité des variances, par traitement et par position du mur :
#     test de Bartlett et rapport entre le plus grand et le plus petit
#     écart-type des résidus (règle usuelle : rapport inférieur à 2).
bartlett.test(residuals(m) ~ dm$ttt)
bartlett.test(residuals(m) ~ dm$mur)
et_ttt <- tapply(residuals(m), dm$ttt, sd)
et_mur <- tapply(residuals(m), dm$mur, sd)
max(et_ttt) / min(et_ttt)
max(et_mur) / min(et_mur)

# 3.3 Additivité : test F de l'interaction traitement x mur
m_int <- lmer(prop_C ~ ttt * mur + (1 | jour), data = dm)
anova(m_int, type = 2)

# 3.4 Normalité des effets aléatoires du jour (sans objet si leur variance
#     est nulle : tous les effets valent alors 0)
ranef(m)

# 3.5 Observations influentes : distance de Cook maximale (seuil usuel : 1),
#     calculée sur le modèle équivalent sans le jour
max(cooks.distance(a))

# Graphiques de diagnostic : résidus vs valeurs ajustées et QQ-plot
par(mfrow = c(1, 2))
plot(fitted(m), residuals(m), xlab = "Valeurs ajustées", ylab = "Résidus")
abline(h = 0, lty = 2)
qqnorm(residuals(m)); qqline(residuals(m))
par(mfrow = c(1, 1))


# 4. TESTS RÉPONDANT AUX QUESTIONS DE RECHERCHE ###############################
# Moyennes ajustées : proportion moyenne de chaque traitement estimée pour la
# situation sans mur (mur = "N"), à partir de l'ensemble des essais.
# Degrés de liberté de Satterthwaite, comme pour le tableau d'ANOVA (section 2).

emm_options(lmer.df = "satterthwaite")
em     <- emmeans(m, ~ ttt, at = list(mur = "N"))
em_alt <- emmeans(m, ~ ttt,
                  at = list(mur = "N", ttt = c("L", "M", "S", "T")))

# 4.1 Préférence : chaque moyenne ajustée comparée à 0.5 (test t).
#     Plantes alternatives : correction de Holm sur 4 tests.
#     Témoin : test unique, sans correction.
res <- as.data.frame(summary(em, infer = c(TRUE, TRUE), null = 0.5))
alt <- res$ttt != "N"
res$p_holm <- NA
res$p_holm[alt] <- p.adjust(res$p.value[alt], method = "holm")
# IC 95 % : corrigés par Šidák pour les 4 plantes alternatives (cohérents
# avec la correction pour tests multiples), non corrigé pour le témoin.
ic_sidak <- as.data.frame(confint(em_alt, adjust = "sidak"))
res$IC_inf <- res$lower.CL
res$IC_sup <- res$upper.CL
res$IC_inf[alt] <- ic_sidak$lower.CL[match(res$ttt[alt], ic_sidak$ttt)]
res$IC_sup[alt] <- ic_sidak$upper.CL[match(res$ttt[alt], ic_sidak$ttt)]
res

# 4.2 Plantes alternatives comparées entre elles (test de Tukey).
#     Une différence négative signifie que la première plante citée attire
#     davantage de Lygus que la seconde (le concombre servant de référence
#     commune).
pairs(em_alt, adjust = "tukey", infer = c(TRUE, TRUE))
cld_alt <- cld(em_alt, Letters = letters, adjust = "tukey")
cld_alt

# 4.3 Chaque plante alternative comparée au témoin (test de Dunnett).
#     p-values exactes ("mvt", calcul numérique de la loi multivariée) ;
#     set.seed rend ce calcul reproductible.
set.seed(2026)
dun <- as.data.frame(summary(contrast(em, method = "trt.vs.ctrl", ref = 1,
                                      adjust = "mvt"),
                             infer = c(TRUE, TRUE)))
dun

# 4.4 Effet du mur : moyennes ajustées (moyennées sur les traitements) et
#     comparaisons deux à deux (Tukey)
em_mur <- emmeans(m, ~ mur)
em_mur
pairs(em_mur, adjust = "tukey", infer = c(TRUE, TRUE))


# 5. TAUX DE RÉPONSE ##########################################################
# Vérifie que la proportion de Lygus ayant choisi une cage ne dépend pas du
# traitement (test de Kruskal-Wallis, distribution asymétrique proche de 1).

tapply(dm$taux_reponse, dm$ttt, mean)
kruskal.test(taux_reponse ~ ttt, data = dm)


# 6. FIGURES ##################################################################

lab_ttt <- c(N = "Témoin (terreau)", L = "Luzerne", M = "Moutarde",
             S = "Sarrasin", T = "Tournesol")
lab_mur <- c(N = "Aucun mur", A = "Mur côté plante\nalternative",
             C = "Mur côté\nconcombre")

theme_pub <- theme_classic(base_size = 9) +
  theme(axis.text  = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        axis.line  = element_line(colour = "black", linewidth = 0.3),
        plot.margin = margin(4, 6, 4, 4),
        legend.position = "none")

# Titre de l'axe des ordonnées, avec Lygus en italique
y_lab <- expression("Proportion de " * italic("Lygus") *
                      " dans la cage du concombre")

# Rappel du sens de lecture sur l'axe de droite
theme_sec <- theme(axis.text.y.right  = element_text(colour = "grey30",
                                                     size = 7.5, hjust = 0),
                   axis.ticks.y.right = element_blank(),
                   axis.line.y.right  = element_blank())

etoiles <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                       ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s."))))

# --- Figure 1 : préférence entre le concombre et les plantes alternatives ---
# (questions 1 à 3 ; le témoin est présenté à part, figure 2)
# Points noirs : moyennes ajustées (sans mur). Barres : IC 95 % corrigés par
# Šidák. Cercles gris : proportion observée dans chaque essai (non corrigée).
# Étoiles : test contre 0.5 (correction de Holm).
# Lettres : test de Tukey entre plantes alternatives ; deux plantes
# partageant une lettre ne diffèrent pas significativement.

alts    <- c("L", "M", "S", "T")
fig_dat <- res[alt, ]
fig_dat$sig <- etoiles(fig_dat$p_holm)
fig_dat <- merge(fig_dat,
                 data.frame(ttt = as.character(cld_alt$ttt),
                            lettre = trimws(cld_alt$.group)),
                 by = "ttt")
fig_dat$ttt <- factor(fig_dat$ttt, levels = alts)
pts_alt <- droplevels(subset(dm, ttt %in% alts))

n_ttt <- table(dm$ttt)
lab_x <- paste0(lab_ttt[alts], "\nn = ", n_ttt[alts])

fig1 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             colour = "grey40") +
  geom_point(data = pts_alt, aes(x = ttt, y = prop_C),
             position = position_jitter(width = 0.15, height = 0, seed = 1),
             shape = 1, size = 1.1, colour = "grey60", stroke = 0.35) +
  geom_errorbar(data = fig_dat,
                aes(x = ttt, ymin = IC_inf, ymax = IC_sup),
                width = 0.12, linewidth = 0.45) +
  geom_point(data = fig_dat, aes(x = ttt, y = emmean),
             shape = 21, fill = "black", size = 2.2) +
  geom_text(data = fig_dat, aes(x = ttt, y = 1.07, label = sig), size = 2.9) +
  geom_text(data = fig_dat, aes(x = ttt, y = 1.15, label = lettre),
            size = 3, fontface = "italic") +
  scale_x_discrete(labels = lab_x) +
  scale_y_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0)),
                     sec.axis = dup_axis(breaks = c(0.8, 0.2), name = NULL,
                       labels = c("Préférence\nconcombre",
                                  "Préférence\nplante opposée"))) +
  labs(x = NULL, y = y_lab) +
  theme_pub + theme_sec
fig1
ggsave("Fig1_preference.pdf", fig1, width = 130, height = 95, units = "mm",
       device = cairo_pdf)
ggsave("Fig1_preference.png", fig1, width = 130, height = 95, units = "mm",
       dpi = 600)

# --- Figure 2 : témoin négatif (questions 4 et 5) ---------------------------
# A : moyenne ajustée du témoin (sans mur) et IC 95 % (test unique) ; cercles
#     gris : essais observés ; étoiles : test contre 0.5.
# B : écart de proportion de Lygus dans la cage de la plante par rapport au
#     terreau, c'est-à-dire l'opposé de la différence « plante - témoin »
#     calculée en 4.3 (une valeur positive signifie que la plante attire plus
#     que le terreau). Barres : IC 95 % simultanés de Dunnett ; étoiles :
#     test de Dunnett.

tem_dat <- res[!alt, ]
tem_dat$sig <- etoiles(tem_dat$p.value)
pts_tem <- droplevels(subset(dm, ttt == "N"))

fig2a <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             colour = "grey40") +
  geom_point(data = pts_tem, aes(x = ttt, y = prop_C),
             position = position_jitter(width = 0.12, height = 0, seed = 1),
             shape = 1, size = 1.1, colour = "grey60", stroke = 0.35) +
  geom_errorbar(data = tem_dat,
                aes(x = ttt, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 0.45) +
  geom_point(data = tem_dat, aes(x = ttt, y = emmean),
             shape = 21, fill = "black", size = 2.2) +
  geom_text(data = tem_dat, aes(x = ttt, y = 1.07, label = sig),
            size = 2.9) +
  scale_x_discrete(labels = paste0("Témoin\n(terreau)\nn = ",
                                   nrow(pts_tem))) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0))) +
  labs(x = NULL, y = y_lab, tag = "A") +
  theme_pub

dun_fig <- data.frame(ttt    = factor(sub(" - N", "", dun$contrast),
                                      levels = alts),
                      ecart  = -dun$estimate,
                      IC_inf = -dun$upper.CL,
                      IC_sup = -dun$lower.CL,
                      sig    = etoiles(dun$p.value))

fig2b <- ggplot(dun_fig, aes(x = ttt, y = ecart)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey40") +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.12,
                linewidth = 0.45) +
  geom_point(shape = 21, fill = "black", size = 2.2) +
  geom_text(aes(y = 0.37, label = sig), size = 2.9) +
  scale_x_discrete(labels = lab_ttt[alts]) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2),
                     sec.axis = dup_axis(breaks = c(0.2, -0.2), name = NULL,
                       labels = c("Plante plus\nattractive que\nle terreau",
                                  "Plante moins\nattractive\nque le terreau"))
                     ) +
  labs(x = NULL, tag = "B",
       y = expression(atop("Écart de proportion de " * italic("Lygus") *
                             " dans la cage",
                           "de la plante par rapport au terreau"))) +
  theme_pub + theme_sec

fig2 <- fig2a + fig2b + plot_layout(widths = c(1, 3.2))
fig2
ggsave("Fig2_temoin.pdf", fig2, width = 160, height = 90, units = "mm",
       device = cairo_pdf)
ggsave("Fig2_temoin.png", fig2, width = 160, height = 90, units = "mm",
       dpi = 600)

# --- Figure S1 : effet du mur ------------------------------------------------
# Points noirs et barres : moyennes ajustées par position du mur et IC 95 %
# corrigés par Šidák (3 intervalles).
# Symboles gris : proportion observée dans chaque essai, par traitement.
# Lettres : test de Tukey entre positions du mur.

cld_mur <- cld(em_mur, Letters = letters, adjust = "tukey")
mur_dat <- merge(as.data.frame(confint(em_mur, adjust = "sidak")),
                 data.frame(mur = cld_mur$mur,
                            lettre = trimws(cld_mur$.group)))
n_mur  <- table(dm$mur)
lab_mx <- paste0(lab_mur, "\nn = ", n_mur[names(lab_mur)])

figS1 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             colour = "grey40") +
  geom_point(data = dm, aes(x = mur, y = prop_C, shape = ttt),
             position = position_jitter(width = 0.15, height = 0, seed = 2),
             size = 1.2, colour = "grey55", stroke = 0.35) +
  geom_errorbar(data = mur_dat,
                aes(x = mur, ymin = lower.CL, ymax = upper.CL),
                width = 0.12, linewidth = 0.45) +
  geom_point(data = mur_dat, aes(x = mur, y = emmean),
             shape = 21, fill = "black", size = 2.2) +
  geom_text(data = mur_dat, aes(x = mur, y = 1.07, label = lettre),
            size = 3, fontface = "italic") +
  scale_shape_manual(values = c(N = 4, L = 1, M = 2, S = 0, T = 5),
                     labels = lab_ttt, name = NULL) +
  scale_x_discrete(labels = lab_mx) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0))) +
  labs(x = NULL, y = y_lab) +
  theme_pub +
  theme(legend.position = "right", legend.key.size = unit(3.5, "mm"))
figS1
ggsave("FigS1_effet_mur.pdf", figS1, width = 140, height = 90, units = "mm",
       device = cairo_pdf)
ggsave("FigS1_effet_mur.png", figS1, width = 140, height = 90, units = "mm",
       dpi = 600)


# 7. ANALYSE COMPLÉMENTAIRE : CONFIGURATION DES CAGES (TOURNESOL) ##########
# Pour le tournesol, deux configurations de cages ont été utilisées (colonne
# grd_cage). Le traitement tournesol est scindé en deux niveaux (T_A, T_AC)
# dans le modèle linéaire mixte, et la différence T_A - T_AC est testée par
# un test t (comparaison unique). Le bruit résiduel et l'effet du mur sont
# estimés sur l'ensemble des essais, ce qui est nécessaire car le mur est en
# partie confondu avec la configuration (les 4 essais tournesol avec le mur
# côté tournesol sont tous en configuration AC).

dm$ttt_cage <- ifelse(dm$ttt == "T", paste0("T_", dm$grd_cage),
                      as.character(dm$ttt))
dm$ttt_cage <- factor(dm$ttt_cage,
                      levels = c("N", "L", "M", "S", "T_A", "T_AC"))
table(Traitement = dm$ttt_cage, Mur = dm$mur)

m_cage  <- lmer(prop_C ~ ttt_cage + mur + (1 | jour), data = dm)
em_cage <- emmeans(m_cage, ~ ttt_cage,
                   at = list(mur = "N", ttt_cage = c("T_A", "T_AC")))

# Moyennes ajustées (sans mur) des deux configurations, comparées à 0.5
summary(em_cage, infer = c(TRUE, TRUE), null = 0.5)

# Différence entre les deux configurations (test t)
diff_cage <- as.data.frame(summary(pairs(em_cage), infer = c(TRUE, TRUE)))
diff_cage

# Vérification sur les seuls essais tournesol (modèle sans le jour)
tt <- droplevels(subset(dm, ttt == "T"))
drop1(lm(prop_C ~ grd_cage + mur, data = tt), test = "F")

# --- Figure S2 : configuration des cages (tournesol) -------------------------
# Points noirs et barres : moyennes ajustées (sans mur) et IC 95 %.
# Symboles gris : proportion observée dans chaque essai, selon la position
# du mur. Crochet : test t de la différence entre configurations.

cage_dat <- as.data.frame(summary(em_cage))
tt$ttt_cage <- factor(paste0("T_", tt$grd_cage), levels = c("T_A", "T_AC"))
n_cage <- table(tt$ttt_cage)
lab_cx <- c(T_A  = "Grande cage :\ntournesol seul",
            T_AC = "Grandes cages :\ntournesol et concombre")
lab_cx <- paste0(lab_cx, "\nn = ", n_cage[names(lab_cx)])
p_txt  <- paste0("p = ", formatC(diff_cage$p.value, format = "f", digits = 2))

figS2 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             colour = "grey40") +
  geom_point(data = tt, aes(x = ttt_cage, y = prop_C, shape = mur),
             position = position_jitter(width = 0.12, height = 0, seed = 3),
             size = 1.3, colour = "grey55", stroke = 0.4) +
  geom_errorbar(data = cage_dat,
                aes(x = ttt_cage, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 0.45) +
  geom_point(data = cage_dat, aes(x = ttt_cage, y = emmean),
             shape = 21, fill = "black", size = 2.2) +
  annotate("segment", x = 1, xend = 2, y = 1.06, yend = 1.06,
           linewidth = 0.3) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 1.06, yend = 1.03,
           linewidth = 0.3) +
  annotate("text", x = 1.5, y = 1.1, label = p_txt, size = 2.9) +
  scale_shape_manual(values = c(N = 1, A = 2, C = 0),
                     labels = c(N = "Aucun mur", A = "Mur côté tournesol",
                                C = "Mur côté concombre"),
                     name = NULL, drop = FALSE) +
  scale_x_discrete(labels = lab_cx) +
  scale_y_continuous(limits = c(0, 1.13), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0))) +
  labs(x = NULL, y = y_lab) +
  theme_pub +
  theme(legend.position = "right", legend.key.size = unit(3.5, "mm"))
figS2
ggsave("FigS2_cages_tournesol.pdf", figS2, width = 130, height = 90,
       units = "mm", device = cairo_pdf)
ggsave("FigS2_cages_tournesol.png", figS2, width = 130, height = 90,
       units = "mm", dpi = 600)


# 8. ENVIRONNEMENT ############################################################
sessionInfo()
