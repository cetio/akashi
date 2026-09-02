module akashi.request;

import core.thread : Thread;
import core.time : dur;
import std.array : join;
import std.datetime : Clock, Duration, SysTime;
import std.net.curl : HTTP;
import std.string : replace;
import std.uri : encode;

package(akashi):

struct RequestClient
{
private:
    string host;
    long minIntervalMs;
    HTTP http;
    bool initialized;
    SysTime lastRequestTime;

    ubyte[] request(
        HTTP.Method method,
        string path,
        string[string] query = null,
        string data = null)
    {
        rateLimit();
        if (!initialized)
        {
            http = HTTP();
            initialized = true;
        }

        ushort status;
        ubyte[] ret;
        http.clearRequestHeaders();
        http.url = buildURL(path, query);
        http.method = method;
        if (method == HTTP.Method.post)
            http.setPostData(data, "application/x-www-form-urlencoded");
        else
            http.onSend = null;

        http.onReceiveStatusLine = (HTTP.StatusLine line) {
            status = line.code;
        };
        http.onReceive = (ubyte[] chunk) {
            if (chunk.length > 0)
                ret ~= chunk;

            return chunk.length;
        };
        http.perform();

        if (status >= 200 && status < 300)
            return ret;
        return null;
    }

    void rateLimit()
    {
        if (lastRequestTime == SysTime.init)
        {
            lastRequestTime = Clock.currTime();
            return;
        }

        Duration minInterval = dur!"msecs"(minIntervalMs);
        Duration elapsed = Clock.currTime() - lastRequestTime;
        if (elapsed < minInterval)
            Thread.sleep(minInterval - elapsed);

        lastRequestTime = Clock.currTime();
    }

    string buildURL(string path, string[string] query = null)
    {
        string ret = host~path;
        if (query.length == 0)
            return ret;

        string[] pairs;
        foreach (key, value; query)
            pairs ~= encode(key).replace("+", "%2B")~"="~encode(value).replace("+", "%2B");

        return ret~"?"~pairs.join("&");
    }

package:
    this(string host, long minIntervalMs)
    {
        this.host = host;
        this.minIntervalMs = minIntervalMs;
    }

    ubyte[] get(string path, string[string] query = null)
    {
        return request(HTTP.Method.get, path, query);
    }

    ubyte[] post(string path, string data)
    {
        return request(
            HTTP.Method.post,
            path,
            null,
            data
        );
    }
}

unittest
{
    RequestClient client = RequestClient("https://example.com", 0);
    assert(client.buildURL("/path") == "https://example.com/path");
    assert(client.buildURL("/path", ["+": "a+b"]) == "https://example.com/path?%2B=a%2Bb");
}
