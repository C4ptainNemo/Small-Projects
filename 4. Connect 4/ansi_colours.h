#ifndef ANSI_COLORS_H
#define ANSI_COLORS_H

// ANSI escape codes for foreground (text) colors
#define FG_BLACK            "\033[30m"  // Black
#define FG_RED              "\033[31m"  // Red
#define FG_GREEN            "\033[32m"  // Green
#define FG_YELLOW           "\033[33m"  // Yellow
#define FG_BLUE             "\033[34m"  // Blue
#define FG_MAGENTA          "\033[35m"  // Magenta
#define FG_CYAN             "\033[36m"  // Cyan
#define FG_WHITE            "\033[37m"  // White

// ANSI escape codes for bright foreground (text) colors
#define FG_BRIGHT_BLACK     "\033[90m"  // Bright Black
#define FG_BRIGHT_RED       "\033[91m"  // Bright Red
#define FG_BRIGHT_GREEN     "\033[92m"  // Bright Green
#define FG_BRIGHT_YELLOW    "\033[93m"  // Bright Yellow
#define FG_BRIGHT_BLUE      "\033[94m"  // Bright Blue
#define FG_BRIGHT_MAGENTA   "\033[95m"  // Bright Magenta
#define FG_BRIGHT_CYAN      "\033[96m"  // Bright Cyan
#define FG_BRIGHT_WHITE     "\033[97m"  // Bright White

// ANSI escape codes for background colors
#define BG_BLACK            "\033[40m"  // Background Black
#define BG_RED              "\033[41m"  // Background Red
#define BG_GREEN            "\033[42m"  // Background Green
#define BG_YELLOW           "\033[43m"  // Background Yellow
#define BG_BLUE             "\033[44m"  // Background Blue
#define BG_MAGENTA          "\033[45m"  // Background Magenta
#define BG_CYAN             "\033[46m"  // Background Cyan
#define BG_WHITE            "\033[47m"  // Background White

// ANSI escape codes for bright background colors
#define BG_BRIGHT_BLACK     "\033[100m" // Bright Background Black
#define BG_BRIGHT_RED       "\033[101m" // Bright Background Red
#define BG_BRIGHT_GREEN     "\033[102m" // Bright Background Green
#define BG_BRIGHT_YELLOW    "\033[103m" // Bright Background Yellow
#define BG_BRIGHT_BLUE      "\033[104m" // Bright Background Blue
#define BG_BRIGHT_MAGENTA   "\033[105m" // Bright Background Magenta
#define BG_BRIGHT_CYAN      "\033[106m" // Bright Background Cyan
#define BG_BRIGHT_WHITE     "\033[107m" // Bright Background White

// ANSI escape codes for text styles
#define RESET               "\033[0m"
#define TXT_BOLD            "\033[1m"   // Bold
#define TXT_DIM             "\033[2m"   // Dim
#define TXT_ITALIC          "\033[3m"   // Italic
#define TXT_UNDERLINE       "\033[4m"   // Underline
#define TXT_BLINK           "\033[5m"   // Blink
#define TXT_INVERSE         "\033[7m"   // Inverse
#define TXT_HIDDEN          "\033[8m"   // Hidden
#define TXT_STRIKETHROUGH   "\033[9m"   // Strikethrough

#endif // ANSI_COLORS_H