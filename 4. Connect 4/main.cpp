#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <array>
#include <cstdlib>

#include "ansi_colours.h"

const size_t boardRows = 6;
const size_t boardColumns = 7;
using gameBoard = std::array<std::array<int, boardColumns>, boardRows>; // Row x Column, values is 0 = empty, 1 = red, 2 = yellow

gameBoard boardState;

void clearConsole()
{
  std::system("cls"); // clears the console
}

bool placePiece(gameBoard &board, size_t column, int playerNumber) // Places a piece in the board. Return true if piece was placed
{
  // Check if the column is within the valid range
  if (column >= board[0].size()) return false;

  // Start from the bottom row
  size_t row = board.size() - 1;

  // Find the first empty row from the bottom up
  while (row < board.size() && board[row][column] != 0) 
    {
        if (row == 0) return false; // If we've reached the top and it's occupied, return false
        row--;
    }

  // Place the piece
  board[row][column] = playerNumber;
  return true;
}

size_t inputColumn(const size_t &maxColumns)
{
  size_t column;
  
  while (1)
  {
    std::cout << "Enter Column Number - 1 to " << maxColumns << ": ";
    if (!(std::cin >> column))
    {
      std::cin.clear();
      std::cin.ignore(10000, '\n');
      std::cerr << "Invalid Input - Enter a number" << std::endl;
      continue;
    }
    if (column > maxColumns)
    {
      std::cout << "Column is out of range..." << std::endl;
      continue;
    }
    return column - 1; // -1 since index is 1 less
  }
}

void drawBoard(gameBoard board)
{
  std::cout << BG_GREEN;
  for (size_t i = 1; i < board[0].size() + 1; i++)
  {
    if (i > 9)
    {
      std::cout << " " << i;
    }
    else
    {
      std::cout << " " << i << " ";
    }
    
  }
  std::cout << RESET << std::endl;
  
  for (auto &&row : board)
  {
    for (auto &&column : row)
    {
      switch (column)
      {
      case 0:
        std::cout << BG_BRIGHT_WHITE << " O " << RESET;
        break;

      case 1:
        std::cout << BG_RED << " O " << RESET;
        break;

      case 2:
        std::cout << BG_YELLOW << " O " << RESET;
        break;
      
      default:
        break;
      }
    }
    std::cout << std::endl;
  }
  
}

bool checkConnect4(gameBoard &board, int playerNumber)
{
  size_t rowMax = board.size();
  size_t columnMax = board[0].size();
  int counter;
  
  // Horizontal Check
  for (size_t row = 0; row < rowMax; row++)
  {
    counter = 0;
    for (size_t column = 0; column < columnMax; column++)
    {
      if (board[row][column] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
    }
  }

  size_t rowCheck;
  size_t columnCheck;
  
  // Vertical Check
  for (size_t column = 0; column < columnMax; column++)
  {
    counter = 0;
    for (size_t row = 0; row < rowMax; row++)
    {
      if (board[row][column] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
    }
  }

  // Diagonal Check (/) Start at Top Left
  for (size_t row = 0; row < rowMax; row++)
  {
    counter = 0;
    rowCheck = row; 
    columnCheck = 0;
    while (1)
    {
      if (board[rowCheck][columnCheck] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
      if (rowCheck == 0 || columnCheck == (columnMax - 1)) break;
      rowCheck--;
      columnCheck++;
    }  
  }
  for (size_t column = 1; column < columnMax; column++)
  {
    counter = 0;
    rowCheck = rowMax; 
    columnCheck = column;
    while (1)
    {
      if (board[rowCheck][columnCheck] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
      if (rowCheck == 0 || columnCheck == (columnMax - 1)) break;
      rowCheck--;
      columnCheck++;
    }
  }

  //  Diagonal Check (\) Start at Bottom Left
  for (size_t row = rowMax - 1; row > 0; row--)
  {
    counter = 0;
    rowCheck = row; 
    columnCheck = 0;
    while (1)
    {
      if (board[rowCheck][columnCheck] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
      if (rowCheck == (rowMax - 1) || columnCheck == (columnMax - 1)) break;
      rowCheck++;
      columnCheck++;
    }  
  }
  for (size_t column = 1; column < columnMax; column++)
  {
    counter = 0;
    rowCheck = 0; 
    columnCheck = column;
    while (1)
    {
      if (board[rowCheck][columnCheck] == playerNumber)
      {
        ++counter;
        if (counter == 4) return true;
      }
      else counter = 0;
      if (rowCheck == (rowMax - 1) || columnCheck == (columnMax - 1)) break;
      rowCheck++;
      columnCheck++;
    }  
  }
  
  
  
  



  return false; // If not connect4
}

int main()
{
  while (1)
  {
    clearConsole();
    drawBoard(boardState);
    std::cout << "Player 1" << std::endl;
    placePiece(boardState, inputColumn(boardColumns), 1);
    if (checkConnect4(boardState, 1))
    {
      std::cout << "Player " << 1 << " Won!" << std::endl;
      break;
    }

    clearConsole();
    drawBoard(boardState);
    std::cout << "Player 2" << std::endl;
    placePiece(boardState, inputColumn(boardColumns), 2);
    if (checkConnect4(boardState, 2))
    {
      std::cout << "Player " << 2 << " Won!" << std::endl;
      break;
    }
  }
  drawBoard(boardState);
  while (1)
  {
    int choice;
    std::cout << "Close Game?\n1. Yes\n2. No\nEnter Number: ";
    std::cin >> choice;
    if (choice == 1) break;
  }
  return 0;
}