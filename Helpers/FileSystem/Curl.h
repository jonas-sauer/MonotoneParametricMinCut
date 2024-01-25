#pragma once

#include <string>
#include <stdio.h>
#include <curl/curl.h>

#include "../Assert.h"

namespace Curl {

    size_t write_data(void* ptr, size_t size, size_t nmemb, FILE *stream) noexcept {
        size_t written = fwrite(ptr, size, nmemb, stream);
        return written;
    }

    void downloadFile(const std::string& url, const std::string fileName) noexcept {
        CURL *curl;
        curl = curl_easy_init();
        if (curl) {
            FILE *fp;
            fp = fopen(fileName.c_str(),"wb");
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
            CURLcode res = curl_easy_perform(curl);
            if (res != 0) {
                error("Error while downloading file!");
            }
            curl_easy_cleanup(curl);
            fclose(fp);
        }
    }

}