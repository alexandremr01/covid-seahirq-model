#include "real_data_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int BUFFER_SIZE = 4096;

// https://stackoverflow.com/questions/12911299/read-csv-file-in-c
const char* getfield(char *line, int num, char separator)
{
    char* first = line;
    char* last = strchr(line, separator);
    if (line[0] == '\"')
    {
        char* first = line + 1;
        last = strchr(first, '\"');
        last++;
    }
    else
    {
        last = strchr(line, separator);
    }

    if (num == 0)
    {
        *last = '\0';
        return line;
    }
    while (num > 0)
    {
        first = last + 1;
        if (*first == '\"')
        {
            last = strchr(first, '\"');
            last++;
        }
        else
        {
            last = strchr(first, separator);
            if (!last) return NULL;
        }
        num--;
    }
    *last = '\0';
    return first;
}

DataReal reader(const char* csv_file, const char* country, char separator)
{
    DataReal result;
    for (int i = 0; i < MAX_DAYS_REAL; ++i)
    {
        result.data[i] = 0.0;
    }
    result.size = 0;
	FILE* csv;
	fopen_s(&csv, csv_file, "r");

    char line[4096];
    char tmp[4096];
    while (fgets(line, 4096, csv) != NULL)
    {
        int i = 4;
        bool finished = false;
        memcpy(tmp, line, 4096);
        const char* token = NULL;
        token = getfield(tmp, 1, separator);
        if (token && strcmp(token, country) == 0)
        {
            while (!finished)
            {
                memcpy(tmp, line, 4096);
                token = getfield(tmp, i, separator);
                if (token == NULL)
                {
                    result.size = i - 4;
                    break;
                }
                result.data[i - 4] += atof(token);
                i++;
            }
        }


        if (feof(csv)) break;
    }
    int i;
    for (i = 0; i < result.size; ++i)
    {
        if (result.data[i] > 0.00001) break;
    }

    memmove(&result.data[0], &result.data[i], (result.size - i) * sizeof(double));
    result.size -= i;

    fclose(csv);

    return result;
}

bool read_non_comment_line(char* buffer, FILE* input) {
  while(fgets(buffer, BUFFER_SIZE, input)) {
    if (buffer[0] != '#' && buffer[0] != '\n')
      return true;
  }
  return false;
}

bool read_simple_array(double array[NEA], char* buffer, FILE* input) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars=0;
  int chars_old;
  double val;
  char title[128];
  sscanf_s(buffer, "%s%n", title, (unsigned)_countof(title), &chars);
  //printf("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf_s(buffer + chars, "%lf%n", &val, &chars_old);
    chars += chars_old;
    array[i] = val;
    //printf("  %lf", val);
  }
  //printf("\n");
  return true;
}

bool read_simple_array(int array[NEA], char* buffer, FILE* input) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars=0;
  int chars_old;
  int val;
  char title[128];
  sscanf_s(buffer, "%s%n", title, (unsigned)_countof(title), &chars);
  //printf("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf_s(buffer + chars, "%d%n", &val, &chars_old);
    chars += chars_old;
    array[i] = val;
    //printf("  %d", val);
  }
  //printf("\n");
  return true;
}

bool read_array_per_day(double array[MAX_DAYS][NEA], char* buffer, FILE* input, int day) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars = 0;
  int chars_old;
  double val;
  char title[128];
  sscanf_s(buffer, "%s%n", title, (unsigned)_countof(title), &chars);
  //printf_s("title: %s\n", title);
  for(int i=0; i < NEA; ++i) {
    sscanf_s(buffer + chars, "%lf%n", &val, &chars_old);
    chars += chars_old;
    for (int k = day; k < MAX_DAYS; ++k)
      {
        int idx = day == 0 ? 0 : k - 1;
        array[k][i] = val  >= 0 ? val  : array[idx][i];
      }
    //printf_s("  %lf", val);
  }
  //printf_s("\n");
  return true;
}

bool read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], char* buffer, FILE* input, int day) {
  if (!read_non_comment_line(buffer, input))
    return false;
  int chars = 0;
  int chars_old;
  double val;
  char title[128];
  sscanf_s(buffer, "%s%n", title, (unsigned)_countof(title), &chars);
  //printf_s("title: %s\n", title);
  for(int j = 0; j < NEA; ++j) {
    if (j>0) {
      read_non_comment_line(buffer, input);
      chars = 0;
    }
    for(int i = 0; i < NEA; ++i) {
      sscanf_s(buffer + chars, "%lf%n", &val, &chars_old);
      //printf_s("  %lf", val);
      chars += chars_old;
      for (int k = day; k < MAX_DAYS; ++k) {
        int idx = day == 0 ? 0 : k - 1;
        matrix[k][j][i] = val  >= 0 ? val  : matrix[idx][j][i];
      }
    }
    //printf_s("\n");
  }
  return true;
}

