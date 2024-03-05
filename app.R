library(shiny)
library(deSolve)
library(tidyverse)
library(ggplot2)

# Define the consumer-resource model

C2R2model <- function(time, state, parms){
  with(as.list(c(state, parms)), {
    dN1_dt = e1*a1a*Ra*N1 + e1*a1b*Rb*N1 - d*N1
    dN2_dt = e2*a2a*Ra*N2 + e2*a2b*Rb*N2 - d*N2
    dRa_dt = d*(Sa-Ra) - (a1a*N1*Ra) - (a2a*N2*Ra) 
    dRb_dt = d*(Sb-Rb) - (a1b*N1*Rb) - (a2b*N2*Rb)
    return(list(c(dN1_dt, dN2_dt, dRa_dt, dRb_dt)))
  })
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Consumer-Resource Model (Two Resources)"),
  sidebarLayout(
    sidebarPanel(
      h4("Parameters"),
      sliderInput("a1a", "Consumption rate (Species 1 eats Resource a)", 0.1, 1, 0.4),
      sliderInput("a1b", "Consumption rate (Species 1 eats Resource b)", 0.1, 1, 0.8),
      sliderInput("a2a", "Consumption rate (Species 2 eats Resource a)", 0.1, 1, 0.6),
      sliderInput("a2b", "Consumption rate (Species 2 eats Resource b)", 0.1, 1, 0.5),
      sliderInput("e1", "Conversion coefficient (Species 1)", 0.5, 2, 1),
      sliderInput("e2", "Conversion coefficient (Species 2)", 0.5, 2, 1),
      sliderInput("d", "Mortality rate", 0.05, 0.2, 0.1),
      sliderInput("Sa", "Resource a supply rate", 0.01, 0.5, 0.3),
      sliderInput("Sb", "Resource b supply rate", 0.01, 0.5, 0.3)
    ),
    mainPanel(
      h4("Formula"),
      verbatimTextOutput("formulas"),
      h4("Plots"),
      plotOutput("ZNGI_plot"),
      plotOutput("time_series_plot")
    )
  )
)

# Shiny server
server <- function(input, output) {
  
  parameters <- reactive({
    list(
      a1a = input$a1a, a1b = input$a1b, a2a = input$a2a, a2b = input$a2b,
      e1 = input$e1, e2 = input$e2,
      d = input$d, 
      Sa = input$Sa, Sb = input$Sb
    )
  })
  
  
  
  output$formulas <- renderText({
      "
      dN1_dt = e1*aa1*Ra*N1 + e1*ab1*Rb*N1 - d*N1
      dN2_dt = e2*aa2*Ra*N2 + e2*ab2*Rb*N2 - d*N2
      dRa_dt = d*(Sa-Ra) - (aa1*N1*Ra) - (aa1*N2*Ra) 
      dRb_dt = d*(Sb-Rb) - (ab1*N1*Rb) - (ab2*N2*Rb)"
  })
  
  # Solve the model and create plots
  output$ZNGI_plot <- renderPlot({
    
    ### Model parameters
    times <- seq(0.1, 2000, by = 0.1)
    state <- c(N1 = 0.05, N2 = 0.05, Ra = parameters()$Sa, Rb = parameters()$Sb)
    
    ### Model application
    pop_size <- ode(func = C2R2model, times = times, y = state, parms = parameters())
    n_row <- nrow(pop_size)
    
    ### Slopes and intercepts of the ZNGI's
    ZNGI_slope_N1 <- -parameters()$a1a / parameters()$a1b
    ZNGI_intercept_N1 <- parameters()$d / (parameters()$e1 * parameters()$a1b)  
    ZNGI_slope_N2 <- -parameters()$a2a / parameters()$a2b
    ZNGI_intercept_N2 <- parameters()$d / (parameters()$e2 * parameters()$a2b)
    
    ### Consumption vectors
    eqilibrium_Ra <- pop_size[n_row, 4]
    eqilibrium_Rb <- pop_size[n_row, 5]
    
    convec_df <- data.frame(x = c(eqilibrium_Ra + 4*parameters()$a1a*eqilibrium_Ra, 
                                  eqilibrium_Ra + 4*parameters()$a2a*eqilibrium_Ra),
                            y = c(eqilibrium_Rb + 4*parameters()$a1b*eqilibrium_Rb, 
                                  eqilibrium_Rb + 4*parameters()$a2b*eqilibrium_Rb),
                            xend = c(eqilibrium_Ra - parameters()$a1a*eqilibrium_Ra, 
                                     eqilibrium_Ra - parameters()$a2a*eqilibrium_Ra),
                            yend = c(eqilibrium_Rb - parameters()$a1b*eqilibrium_Rb, 
                                     eqilibrium_Rb - parameters()$a2b*eqilibrium_Rb),
                            species = c("N1", "N2"))
    
    eqilibrium_1 <- pop_size[n_row, 2] %>% as.numeric
    eqilibrium_2 <- pop_size[n_row, 3] %>% as.numeric
    if(eqilibrium_1 < 1e-6){
      convec_df$x[1] <- eqilibrium_Ra
      convec_df$xend[1] <- eqilibrium_Ra
      convec_df$y[1] <- eqilibrium_Rb
      convec_df$yend[1] <- eqilibrium_Rb
    } 
    if(eqilibrium_2 < 1e-6){
      convec_df$x[2] <- eqilibrium_Ra
      convec_df$xend[2] <- eqilibrium_Ra
      convec_df$y[2] <- eqilibrium_Rb
      convec_df$yend[2] <- eqilibrium_Rb
    }
    
    ### Phase diagram
      ggplot() + 
        geom_abline(slope = ZNGI_slope_N1, intercept = ZNGI_intercept_N1, color = "#377EB8", size = 1.2) + 
        geom_abline(slope = ZNGI_slope_N2, intercept = ZNGI_intercept_N2, color = "#E41A1C", size = 1.2) + 
        geom_segment(data = convec_df, aes(x = x, y = y, xend = xend, yend = yend, color = species), linetype = "blank") + 
        geom_segment(data = convec_df, aes(x = x, y = y, xend = xend, yend = yend, color = species), size = 0.5, linetype = "dashed", arrow = arrow(type = "closed", length = unit(0.1, "inches")), show.legend = F) +
        geom_path(data = as.data.frame(pop_size), aes(x = Ra, y = Rb), size = 1.2) +
        geom_point(data = as.data.frame(pop_size), aes(x = last(Ra), y = last(Rb)), size = 2.5) +
        theme_classic(base_size = 14) +
        labs(x = expression(italic(R[a])), y = expression(italic(R[b]))) +
        scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
        scale_color_brewer(name = NULL, palette = "Set1", direction = -1,
                           guide = guide_legend(override.aes = list(
                             linetype = "solid", size = 1.2))) +
        coord_fixed(ratio = 1)
    
  })
  
  output$time_series_plot <- renderPlot({### Model specification
    
    ### Model parameters
    times <- seq(0.1, 2000, by = 0.1)
    state <- c(N1 = 0.05, N2 = 0.05, Ra = parameters()$Sa, Rb = parameters()$Sb)
    
    ### Model application
    pop_size <- ode(func = C2R2model, times = times, y = state, parms = parameters())
    n_row <- nrow(pop_size)
    
    ### Visualize the population dynamics
    pop_size %>%
      as.data.frame() %>%
      gather(key = "Species", value = "N", N1:Rb) %>%
      mutate(trophic = case_when(Species %in% c("N1", "N2") ~ "Consumer",
                                 TRUE ~ "Resource")) %>%
      ggplot(aes(x = time, y = N, color = Species)) + 
      geom_line(size = 1.5) +
      facet_wrap(~ trophic, 
                 ncol = 2, 
                 scales = "free_y",
                 strip.position = "left") +
      theme_classic(base_size = 14) +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "top",
            legend.title = element_blank(),
            plot.margin = margin(r = 8)) + 
      labs(x = "Time", y = NULL) +
      scale_x_continuous(limits = c(0, 2050), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
      scale_color_manual(name = NULL, values = c("blue", "red", "green", "purple"))
    
    
  })
}


# Run the Shiny app
shinyApp(ui, server)

#
# setwd("C:\\Users\\andyh\\OneDrive\\Documents\\shiny\\consumer-resource_two-resources")
# rsconnect::setAccountInfo(name='katasuke', token='318EFBE35FEC40CB6ED4DB2DF70DFBB6', secret='KHl0Bz6aAStTWm3B6cVaC6bnwJiFOAcjsUZ5dg9q')
#rsconnect::deployApp()
