#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

static void ui_draw_extras_limit_speed(UIState *s)
{
    const UIScene *scene = &s->scene;
    cereal::CarControl::SccSmoother::Reader scc_smoother = scene->car_control.getSccSmoother();
    int limit_speed = scc_smoother.getRoadLimitSpeed();
    int left_dist = scc_smoother.getRoadLimitSpeedLeftDist();

    if(limit_speed > 10 && left_dist > 0)
    {
        int w = 200;
        int h = 200;
        int x = (s->scene.viz_rect.x + (bdr_s*2)) + 300;
        int y = 80;
        char str[32];

        nvgBeginPath(s->vg);
        nvgRoundedRect(s->vg, x, y, w, h, 210);
        nvgStrokeColor(s->vg, nvgRGBA(255, 0, 0, 200));
        nvgStrokeWidth(s->vg, 30);
        nvgStroke(s->vg);

        NVGcolor fillColor = nvgRGBA(0, 0, 0, 50);
        nvgFillColor(s->vg, fillColor);
        nvgFill(s->vg);

        nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 250));

        nvgFontSize(s->vg, 140);
        nvgFontFaceId(s->vg, s->font_sans_bold);
        nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

        snprintf(str, sizeof(str), "%d", limit_speed);
        nvgText(s->vg, x+w/2, y+h/2, str, NULL);

        nvgFontSize(s->vg, 120);

        if(left_dist >= 1000)
            snprintf(str, sizeof(str), "%.1fkm", left_dist / 1000.f);
        else
            snprintf(str, sizeof(str), "%dm", left_dist);

        nvgText(s->vg, x+w/2, y+h + 70, str, NULL);
    }
}

static void ui_draw_extras(UIState *s)
{
    ui_draw_extras_limit_speed(s);
}