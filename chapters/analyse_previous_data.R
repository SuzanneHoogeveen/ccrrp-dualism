# Simple script to get the state and context effects from previous studies, 
# including different cultural setting. These will be used to inform the prior
# distributions for our data. 

get_previous_data <- function(){
  # study 1: Watson-Jones 
  va_m_s <- 1.15
  va_b_s <- 0.28
  va_m_t <- 2.98
  va_b_t <- 4.03
  us_m_s <- 1.78
  us_b_s <- 0.74
  us_m_t <- 3.81
  us_b_t <- 2.40
  
  s1 <- data.frame(means = c(va_m_s, va_b_s, va_m_t, va_b_t, 
                             us_m_s, us_b_s, us_m_t, us_b_t),
                   site = rep(c("Vanuatu", "USA"), each=4),
                   state = rep(c("mind","body"), 4),
                   context = rep(c("secular","secular","theistic","theistic"),2))
  s1$prob <- s1$means/7
  
  s1_state <- s1 %>%
    group_by(site,context) %>%
    dplyr::summarise(state_eff = -diff(prob), 
                     p_cont = mean(prob)) %>%
    ungroup(context) %>%
    dplyr::summarise(state_diff = mean(state_eff),
                     p_cont = mean(p_cont))
  
  s1_con <- s1 %>%
    group_by(site,state) %>%
    dplyr::summarise(context_eff = diff(prob),
                     p_cont = mean(prob)) %>%
    ungroup(state) %>%
    dplyr::summarise(context_diff = mean(context_eff),
                     p_cont = mean(p_cont))
  
  s1 <- data.frame(site=s1_state$site, 
                   state_diff = s1_state$state_diff, 
                   context_diff = s1_con$context_diff,
                   p_continues = s1_state$p_cont)
  
  # study 2: Astuti & Harris (2008)
  ma_m <- (7-4.80)/7
  ma_b <- (7-6.22)/7
  
  ma_s <- (14-12.11)/14
  ma_t <- (14-9.92)/14
  
  s2 <- data.frame(site="Madagaskar",state_diff =ma_m-ma_b, 
                   context_diff = ma_t-ma_s,
                   p_continues = mean(c(ma_m,ma_b)))
  
  # study 3: Bering (2002)
  ca_b <- 1-mean(c(.949,.942,.725))
  ca_m <- 1-mean(c(.307,.301,.263))
  
  s3 <- data.frame(site="Florida", 
                   state_diff = ca_m-ca_b, 
                   context_diff = NA, 
                   p_continues = mean(c(ca_m,ca_b)))
  
  # study 4: Barrett et al. (2021)
  s4_dat <- read.csv('../data/BarrettEtAl2021_data.csv')
  
  s4 <- s4_dat %>%
    filter(SLPDEATH == "Death") %>%
    group_by(SITENAME,PRIME,PSYCHPHYS) %>%
    dplyr::summarise(prob = mean(RESPONSE, na.rm=TRUE))
  
  s4_state <- s4 %>%
    group_by(SITENAME,PRIME) %>%
    dplyr::summarise(state_eff = diff(prob), 
                     p_cont = mean(prob)) %>%
    ungroup(PRIME) %>%
    dplyr::summarise(state_diff = mean(state_eff),
                     p_cont = mean(p_cont))
  
  s4_con <- s4 %>%
    group_by(SITENAME,PSYCHPHYS) %>%
    dplyr::summarise(context_eff = diff(prob)) %>%
    ungroup(PSYCHPHYS) %>%
    dplyr::summarise(context_diff = mean(context_eff))
  
  s4 <- data.frame(site=s4_state$SITENAME, 
                   state_diff = s4_state$state_diff, 
                   context_diff = c(NA,s4_con$context_diff),
                   p_continues = s4_state$p_cont)
  
  prevdat <- rbind(s1,s2,s3,s4)
  return(prevdat)
}



