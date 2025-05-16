# Shiny app supplement for Sommer et al 

# Libraries
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(here)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

# Load data
safe_load_data <- function(file_path) {
  tryCatch({
    full_path <- here::here(file_path)
    if (file.exists(full_path)) {
      return(read_csv(full_path))
    }
        if (file.exists(file_path)) {
      return(read_csv(file_path))
    }
        app_dir <- getwd()
    rel_path <- file.path(app_dir, file_path)
    if (file.exists(rel_path)) {
      return(read_csv(rel_path))
    }
    
    stop(paste("File not found:", file_path))
  }, error = function(e) {
    message("Error loading data: ", e$message)
    return(NULL)
  })
}

# Load RDS
safe_load_rds <- function(file_path) {
  tryCatch({
    full_path <- here::here(file_path)
    if (file.exists(full_path)) {
      return(readRDS(full_path))
    }
        if (file.exists(file_path)) {
      return(readRDS(file_path))
    }
        app_dir <- getwd()
    rel_path <- file.path(app_dir, file_path)
    if (file.exists(rel_path)) {
      return(readRDS(rel_path))
    }
    
    stop(paste("File not found:", file_path))
  }, error = function(e) {
    message("Error loading RDS: ", e$message)
    return(NULL)
  })
}

shiny_data <- safe_load_data("shiny_data.csv")
site_comparison_data <- safe_load_data("site_comparison_data.csv")
variable_list <- safe_load_rds("variable_list.rds")

if (is.null(shiny_data) || is.null(site_comparison_data) || is.null(variable_list)) {
  stop("Failed to load required data files. Please ensure all data files are present in the ShinyApp directory.")
}

# Mapping between display names and actual column names
variable_mapping <- list()
for (category in names(variable_list)) {
  for (var_name in names(variable_list[[category]])) {
    display_name <- variable_list[[category]][[var_name]]
    variable_mapping[[display_name]] <- var_name
  }
}

# Define units
variable_units <- list(
  "SORU_Biomass" = "Goldenrod Biomass (g)",
  "POPRC_Biomass" = "Grass Biomass (g)",
  "PlantDiversity" = "Plant Diversity (Shannon-Wiener)",
  "PercentC_SOIL" = "Soil %C",
  "PercentN_SOIL" = "Soil %N",
  "CO2CperHourperg" = "SIR (CO2 per hour per g)",
  "Overall_mineralization_rate" = "N-min (mg N/cm² per month)",
  "PercentC_LITTER" = "Litter %C",
  "PercentN_LITTER" = "Litter %N",
  "PercentC_SORU" = "Goldenrod %C",
  "PercentN_SORU" = "Goldenrod %N",
  "PercentC_POPRC" = "Grass %C",
  "PercentN_POPRC" = "Grass %N"
)

treatment_colors <- c(
  "Reactor" = "#663366",    # Muted Purple
  "Resistor" = "#D4AF37",  # Muted Gold
  "Vegetation" = "#9BA48C"  # Green
)

# UI definition
ui <- fluidPage(
  div(style = "text-align: center; margin-bottom: 30px;",
      h1("Herbivore Populations Restructure Ecosystems", 
         style = "font-size: 28px; margin-bottom: 10px;"),
      h3(HTML("Data from Sommer et al. Herbivore population differences rival geographic and biophysical variation in structuring ecosystem function. <i>Global Change Biology</i>"),
         style = "font-size: 16px; font-weight: normal;")
  ),
  
  div(style = "display: flex; justify-content: center; margin-bottom: 20px;",
      div(style = "width: 50%;",
          selectInput("section", "Select Visualization:",
                     choices = c("Raw Data Visualization", "Common Garden Comparison"),
                     width = "100%")
      )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # Graph controls
      uiOutput("section_controls")
    ),
    
    mainPanel(
      width = 9,
      uiOutput("main_content")
    )
  )
)

# Server
server <- function(input, output, session) {
  output$section_controls <- renderUI({
    switch(input$section,
           "Common Garden Comparison" = {
             tagList(
               selectInput("metric", "Select Metric:",
                          choices = variable_mapping,
                          selected = variable_mapping[[1]]),
               checkboxInput("show_anova", "Show ANOVA (model effects)", value = FALSE),
               checkboxInput("show_summary", "Show summary", value = FALSE)
             )
           },
           "Raw Data Visualization" = {
             tagList(
               selectInput("x_var", "X Variable:",
                          choices = variable_mapping,
                          selected = "SORU_Biomass"),
               selectInput("y_var", "Y Variable:",
                          choices = variable_mapping,
                          selected = "POPRC_Biomass"),
               checkboxGroupInput("treatment_filter", "Select Treatments:",
                                choices = c("Reactor", "Resistor", "Vegetation"),
                                selected = c("Reactor", "Resistor", "Vegetation")),
               radioButtons("year_filter", "Select Year:",
                          choices = c("2021", "2023", "Both"),
                          selected = "Both"),
               selectInput("site_filter", "Select Sites:",
                          choices = unique(shiny_data$Site),
                          selected = unique(shiny_data$Site),
                          multiple = TRUE),
               checkboxInput("show_regression", "Show regression line with confidence interval", value = TRUE)
             )
           }
    )
  })
  
  output$main_content <- renderUI({
    switch(input$section,
           "Common Garden Comparison" = {
             out <- list(
               plotOutput("site_comparison_plot", height = "600px")
             )
             if (!is.null(input$show_anova) && input$show_anova) {
               out <- c(out,
                        list(
                          br(),
                          h4("Does the effect of treatment differ by site and year?"),
                          div("ANOVA table calculated from linear mixed effects models (lmer) with the structure: Metric ~ Site * Treatment * Year + (1|Sample_ID). A non-significant three-way interaction indicates that the effect of treatment is the same across sites and years.", style = "font-style: italic; margin-bottom: 8px;"),
                          DTOutput("anova_table")
                        ))
             }
             if (!is.null(input$show_summary) && input$show_summary) {
               out <- c(out,
                        list(
                          br(),
                          h4("Summary Table"),
                          DTOutput("stats_table")
                        ))
             }
             do.call(tagList, out)
           },
           "Raw Data Visualization" = {
             plotOutput("raw_data_plot", height = "600px")
           }
    )
  })
  
  # stats functions
  perform_stats <- function(data, selected_var) {
    tryCatch({
      # Fit mixed model with three-way interaction
      model <- lmer(as.formula(paste(selected_var, "~ Site * Treatment_PopType * Year + (1|Sample_ID)")),
                   data = data)
      
      emm <- emmeans(model, ~ Treatment_PopType | Site * Year)
      pairs <- pairs(emm, adjust = "tukey")
      
      model_summary <- summary(model)
      
      anova_table <- anova(model)
      
      return(list(
        model = model,
        pairs = pairs,
        summary = model_summary,
        anova = anova_table
      ))
    }, error = function(e) {
      message("Error in statistical analysis: ", e$message)
      return(NULL)
    })
  }
  
  # Site comparison 
  output$site_comparison_plot <- renderPlot({
    req(input$metric)
    
    selected_var <- input$metric
    display_name <- variable_units[[selected_var]]
    
    # summary statistics
    summary_data <- site_comparison_data %>%
      group_by(Site, Treatment_PopType, Year) %>%
      summarise(
        mean = mean(.data[[selected_var]], na.rm = TRUE),
        se = sd(.data[[selected_var]], na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      mutate(Year = case_when(
        Year == "2021" ~ "Initial (2021)",
        Year == "2023" ~ "Final (2023)"
      ))
    
    p <- ggplot() +
      geom_jitter(data = site_comparison_data %>%
                   mutate(Year = case_when(
                     Year == "2021" ~ "Initial (2021)",
                     Year == "2023" ~ "Final (2023)"
                   )),
                 aes(x = Year, y = .data[[selected_var]], 
                     color = Treatment_PopType),
                 width = 0.1, height = 0, alpha = 0.5, size = 3) +
      geom_line(data = summary_data,
                aes(x = Year, y = mean, 
                    color = Treatment_PopType, group = Treatment_PopType),
                size = 1.5) +
      geom_errorbar(data = summary_data,
                   aes(x = Year, y = mean,
                       ymin = mean - se, ymax = mean + se,
                       color = Treatment_PopType),
                   width = 0.1, size = 1) +
      scale_color_manual(values = treatment_colors) +
      labs(x = "Year", 
           y = display_name,
           color = "Treatment") +
      facet_wrap(~Site) +
      scale_x_discrete(limits = c("Initial (2021)", "Final (2023)")) +
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        plot.title = element_text(face = "bold", size = 20)
      )
    
    return(p)
  })
  
  # stat sumary table
  output$stats_table <- renderDT({
    req(input$metric, input$show_summary)
    if (input$show_summary) {
      selected_var <- input$metric
      summary_stats <- site_comparison_data %>%
        group_by(Site, Treatment_PopType, Year) %>%
        summarise(
          mean = mean(.data[[selected_var]], na.rm = TRUE),
          se = sd(.data[[selected_var]], na.rm = TRUE) / sqrt(n()),
          n = n(),
          .groups = "drop"
        ) %>%
        mutate(
          Year = case_when(
            Year == "2021" ~ "Initial (2021)",
            Year == "2023" ~ "Final (2023)"
          ),
          mean = sprintf("%.3f", mean),
          se = sprintf("%.3f", se)
        ) %>%
        arrange(Site, Treatment_PopType, Year) %>%
        transmute(
          Site = Site,
          Treatment = Treatment_PopType,
          Year = Year,
          Mean = mean,
          SE = se,
          n = as.character(n)
        )
      datatable(summary_stats,
                options = list(pageLength = 12,
                             dom = 't',
                             ordering = FALSE),
                rownames = FALSE)
    }
  })
  
  output$anova_table <- renderDT({
    req(input$metric, input$show_anova)
    if (input$show_anova) {
      selected_var <- input$metric
      stats_results <- perform_stats(site_comparison_data, selected_var)
      if (!is.null(stats_results)) {
        anova_table <- as.data.frame(stats_results$anova)
        anova_table$Term <- rownames(anova_table)
        num_cols <- sapply(anova_table, is.numeric)
        anova_table[num_cols] <- lapply(anova_table[num_cols], function(x) sprintf("%.3f", x))
        highlight_row <- grep("Site:Treatment_PopType:Year", anova_table$Term)
        if (length(highlight_row) > 0) {
          anova_table <- rbind(anova_table[highlight_row, ], anova_table[-highlight_row, ])
          highlight_row <- 1
        }
        datatable(anova_table[, c("Term", setdiff(names(anova_table), "Term"))],
                  options = list(pageLength = 10, dom = 't', ordering = FALSE,
                                 rowCallback = JS(
                                   sprintf('function(row, data, index) { if (index == %d) { $(row).css("font-weight", "bold"); $(row).css("background", "#ffffcc"); } }',
                                           highlight_row - 1)
                                 )
                ),
                  rownames = FALSE)
      }
    }
  })
  
  # Raw data
  output$raw_data_plot <- renderPlot({
    req(input$x_var, input$y_var)
    
    x_var <- input$x_var
    y_var <- input$y_var
    x_display <- variable_units[[x_var]]
    y_display <- variable_units[[y_var]]
    
    plot_data <- shiny_data %>%
      filter(Treatment_PopType %in% input$treatment_filter,
             Site %in% input$site_filter,
             Year %in% if(input$year_filter == "Both") c("2021", "2023") else input$year_filter) %>%
      mutate(x = .data[[x_var]], 
             y = .data[[y_var]],
             Year = factor(Year))  # Convert Year to factor
    
    p <- ggplot(plot_data, aes(x = x, y = y, color = Treatment_PopType, shape = Year)) +
      geom_point(alpha = 0.7, size = 3, stroke = 1) +
      scale_color_manual(values = treatment_colors) +
      scale_shape_manual(values = c("2021" = 16, "2023" = 17)) +
      labs(x = x_display,
           y = y_display,
           color = "Treatment",
           shape = "Year",
           title = paste("Relationship between", x_display, "and", y_display),
           subtitle = if(input$year_filter == "Both") {
             "Data from both 2021 and 2023"
           } else {
             paste("Data from", input$year_filter)
           }) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        legend.position = "top",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box = "horizontal",
        legend.spacing.x = unit(1, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        # Margins
        plot.margin = margin(20, 20, 20, 20)
      )
    
    if (input$show_regression) {
      if (input$year_filter == "Both") {
        p <- p + geom_smooth(method = "lm", se = TRUE, 
                            aes(group = Treatment_PopType),
                            alpha = 0.2, size = 1,
                            show.legend = FALSE)
      } else {
        p <- p + geom_smooth(method = "lm", se = TRUE, 
                            aes(group = interaction(Treatment_PopType, Year)),
                            alpha = 0.2, size = 1,
                            show.legend = FALSE)
      }
    }
    
    return(p)
  })
}

# Run the app
shinyApp(ui = ui, server = server) 
